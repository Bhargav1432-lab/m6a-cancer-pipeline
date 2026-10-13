import requests
import json
import statistics

gene_entrez = {
    'CD24': 100133941,
    'AHCY': 191,
    'UBE3C': 9690,
    'RPS6KA3': 6197,
    'CTCF': 10664,
    'FADD': 8772
}

studies = {
    'luad_tcga_pan_can_atlas': 'Lung Adenocarcinoma (TCGA LUAD)',
    'lihc_tcga_pan_can_atlas': 'Liver HCC (TCGA LIHC)'
}

base = 'https://www.cbioportal.org/api'
json_hdr = {'accept': 'application/json', 'Content-Type': 'application/json'}
get_hdr  = {'accept': 'application/json'}

final_results = {}

for study_id, study_name in studies.items():
    print(f"\n=== {study_name} ===")

    # Get molecular profiles — returns a list of profile objects (dicts)
    r = requests.get(f'{base}/studies/{study_id}/molecular-profiles',
                     headers=get_hdr, timeout=30)
    profiles = r.json()   # list of dicts
    print(f"  Available profiles ({len(profiles)}):")
    for p in profiles:
        print(f"    {p['molecularProfileId']}")

    # Pick best expression profile
    priority = ['rna_seq_v2_mrna_median_zscores',
                'mrna_median_zscores',
                'rna_seq_mrna_median_zscores',
                'rna_seq_mrna']
    expr_profile = None
    for kw in priority:
        for p in profiles:
            if kw in p['molecularProfileId']:
                expr_profile = p['molecularProfileId']
                break
        if expr_profile:
            break
    if not expr_profile:
        for p in profiles:
            if 'mrna' in p['molecularProfileId']:
                expr_profile = p['molecularProfileId']
                break

    print(f"  Selected profile: {expr_profile}")

    # Get sample lists
    r2 = requests.get(f'{base}/studies/{study_id}/sample-lists',
                      headers=get_hdr, timeout=30)
    sample_lists = r2.json()
    all_samples_id = None
    for sl in sample_lists:
        if sl['sampleListId'].endswith('_all'):
            all_samples_id = sl['sampleListId']
            break
    print(f"  Sample list: {all_samples_id}")

    if not (expr_profile and all_samples_id):
        print("  Skipping — missing profile or sample list")
        continue

    # Fetch molecular data
    payload = {
        "entrezGeneIds": list(gene_entrez.values()),
        "sampleListId": all_samples_id
    }
    url = f"{base}/molecular-profiles/{expr_profile}/molecular-data/fetch"
    r3 = requests.post(url, headers=json_hdr,
                       data=json.dumps(payload), timeout=90)

    if r3.status_code != 200:
        print(f"  HTTP {r3.status_code}: {r3.text[:200]}")
        continue

    data = r3.json()
    print(f"  Records: {len(data)}")

    # Group by entrez ID -> gene symbol
    gene_vals = {g: [] for g in gene_entrez}
    for rec in data:
        eid = rec.get('entrezGeneId')
        val = rec.get('value')
        sym = next((g for g, e in gene_entrez.items() if e == eid), None)
        if sym and val is not None:
            try:
                gene_vals[sym].append(float(val))
            except (ValueError, TypeError):
                pass

    gene_summary = {}
    print(f"\n  Gene expression z-scores:")
    for gene, vals in gene_vals.items():
        if vals:
            med = statistics.median(vals)
            mn  = min(vals)
            mx  = max(vals)
            n   = len(vals)
            pct_up = sum(1 for v in vals if v > 0) / n * 100
            print(f"    {gene:10s}: n={n:4d}  median z={med:+.3f}  [{mn:.2f}, {mx:.2f}]  {pct_up:.0f}% above median")
            gene_summary[gene] = {
                'n': n,
                'median_zscore': round(med, 3),
                'min': round(mn, 3),
                'max': round(mx, 3),
                'pct_upregulated': round(pct_up, 1)
            }
        else:
            print(f"    {gene:10s}: No data")
            gene_summary[gene] = None

    final_results[study_id] = {
        'study_name': study_name,
        'expr_profile': expr_profile,
        'n_patients': len(set(rec.get('sampleId','') for rec in data)),
        'gene_expression': gene_summary
    }

# Save
with open('tcga_expression_results.json', 'w') as f:
    json.dump(final_results, f, indent=2)

# Print comparison table
print("\n\n=== FINAL COMPARISON TABLE ===")
print(f"{'Gene':<12} {'LUAD median z':>15} {'LIHC median z':>15}  Direction")
print("-" * 65)
genes = list(gene_entrez.keys())
for gene in genes:
    luad_d = final_results.get('luad_tcga_pan_can_atlas', {}).get('gene_expression', {}).get(gene)
    lihc_d = final_results.get('lihc_tcga_pan_can_atlas', {}).get('gene_expression', {}).get(gene)
    lz = f"{luad_d['median_zscore']:+.3f}" if luad_d else "N/A"
    hz = f"{lihc_d['median_zscore']:+.3f}" if lihc_d else "N/A"
    direction = ""
    if luad_d and lihc_d:
        diff = luad_d['median_zscore'] - lihc_d['median_zscore']
        direction = f"Lung {'HIGHER' if diff > 0 else 'LOWER'} by {abs(diff):.3f} z"
    print(f"  {gene:<10} {lz:>15} {hz:>15}  {direction}")

print("\nSaved -> tcga_expression_results.json")
