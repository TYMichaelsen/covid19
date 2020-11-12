def create_location_division_mapping(df, path):
    division_location_mappings = df.groupby(['division', 'location']).groups.keys()
    dl_df = pd.DataFrame(columns=['division','location'])

    for (d,l) in division_location_mappings:
        if d == '' or l == '' or type(d) != str or type(l) != str:
            continue
        dl_df = dl_df.append({
            'division':d,
            'location':l
        }, ignore_index=True)
    dl_df.to_csv(path, sep='\t', index=False)