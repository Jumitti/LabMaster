from pages.design_primer_API import NCBIdna

print(NCBIdna('NM_003130.4', all_slice_forms=False).find_sequences())
print(NCBIdna('NM_003130.4', all_slice_forms=True).find_sequences())