import pandas as pd
import pickle
from tqdm import tqdm


def prepare_search_results(result_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_best_raw.tsv',
                           out_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_web.parquet'):

    print(f"Preparing search results from {result_path} to {out_path}")
    alldf = pd.read_csv(result_path, sep='\t', low_memory=False)
    
    alldf['ref_1_prec_int'] = alldf.apply(lambda x: max(x['ref_1_prec_frag_int'], x['ref_1_prec_frag_water_loss_int']), axis=1)
    
    # drop cols ref_1_prec_frag_int, ref_1_prec_frag_water_loss_int, ref_1_prec_nl_int, ref_1_prec_nl_water_loss_int, total_usage
    alldf = alldf.drop(columns=['ref_1_prec_frag_int', 'ref_1_prec_frag_water_loss_int', 'ref_1_prec_nl_int',
                               'ref_1_prec_nl_water_loss_int', 'total_usage']).reset_index(drop=True)
    
    alldf = alldf[(alldf['ref_1_score'] >= 0.8) & (alldf['ref_1_peak'] >= 3) & 
                  (alldf['ref_1_usage'] >= 0.30) & (alldf['ref_1_prec_int'] >= 0.10)].reset_index(drop=True)
    
    df_1 = alldf[(alldf['ref_2_score'] >= 0.6) & (alldf['ref_2_peak'] >= 3)].reset_index(drop=True)
    df_1['spec_spec'] = True
    
    df_2 = alldf[(alldf['ref_2_id'].isna())].reset_index(drop=True)
    df_2['spec_spec'] = False
    
    df = pd.concat([df_1, df_2], axis=0).reset_index(drop=True)
    del alldf
    
    # load metadata
    print("Loading metadata...")
    # metadata = pd.read_pickle('/Users/shipei/Documents/projects/conjugated_metabolome/db/ms2db/all/all_ms2db_metadata.pkl')
    metadata = pd.read_pickle('/home/shipei/projects/revcos/search/all_ms2db_metadata.pkl')
    db_id_to_db = metadata.set_index('db_id')['db'].to_dict()
    db_id_to_mass = metadata.set_index('db_id')['monoisotopic_mass'].to_dict()
    db_id_to_inchi14 = metadata.set_index('db_id')['inchikey_14'].to_dict()
    
    # add metadata to df
    df['ref_1_mass'] = df['ref_1_id'].map(db_id_to_mass)
    df['delta_mass'] = df['qry_mz'] - df['ref_1_mass'] - 1.007276466812
    df['ref_1_db'] = df['ref_1_id'].map(db_id_to_db)
    df['ref_1_db'] = df['ref_1_db'].apply(lambda x: 1 if x == 'gnps' or x == 'mona' else 0)
    
    def _fill_ref_2_inchi(row):
        if pd.isnull(row['ref_2_id']):
            return str(round(row['delta_mass'], 2))
        else:
            id = str(row['ref_2_id'])
            if id.endswith('.0'):
                id = id[:-2]
            return db_id_to_inchi14[id]
    def _fill_ref_2_db(row):
        if pd.isnull(row['ref_2_id']):
            return 0
        else:
            id = str(row['ref_2_id'])
            if id.endswith('.0'):
                id = id[:-2]
            db = db_id_to_db[id]
            if db == 'gnps' or db == 'mona':
                return 1
            else:
                return 0
    df['ref_2_inchi'] = df.apply(_fill_ref_2_inchi, axis=1)
    df['ref_2_db'] = df.apply(_fill_ref_2_db, axis=1)
    
    df = df[df['delta_mass'] > 1.5].reset_index(drop=True)
    
    # sort by ref_1_inchi, ref_2_inchi, ref_1_db, ref_2_db, ref_1_prec_int, ref_1_usage, ref_1_score, ref_2_score
    df = df.sort_values(by=['ref_1_inchi', 'ref_2_inchi', 'ref_1_db', 'ref_2_db', 'ref_1_prec_int', 'ref_1_usage',
                            'ref_1_score', 'ref_2_score'], ascending=[True, True, False, False, False, False, False, False])
    df = df.reset_index(drop=True)
    
    # def remove_prefix_ref_1(row):
    #     if pd.isnull(row['ref_1_db']):
    #         return None
    #     if row['ref_1_db'] == 'gnps':
    #         return row['ref_1_id'].replace('CCMSLIB000', '')
    #     elif row['ref_1_db'] == 'massbank':
    #         return row['ref_1_id'].replace('MSBNK-', '')
    #     else:
    #         return row['ref_1_id']
    # df['ref_1_id'] = df.apply(remove_prefix_ref_1, axis=1)
    
    # def remove_prefix_ref_2(row):
    #     if pd.isnull(row['ref_2_db']):
    #         return None
    #     if row['ref_2_db'] == 'gnps':
    #         return row['ref_2_id'].replace('CCMSLIB000', '')
    #     elif row['ref_2_db'] == 'massbank':
    #         return row['ref_2_id'].replace('MSBNK-', '')
    #     else:
    #         return row['ref_2_id']
    # df['ref_2_id'] = df.apply(remove_prefix_ref_2, axis=1)
    
    # First get the group counts
    print("Grouping by ref_1_inchi and ref_2_inchi...")
    group_counts = df.groupby(['ref_1_inchi', 'ref_2_inchi']).size().reset_index(name='count')

    # Aggregation
    agg_df = df.groupby(['ref_1_inchi', 'ref_2_inchi']).agg({
        'qry_id': 'first',
        'qry_mz': 'first',
        'ref_1_id': 'first',
        'ref_2_id': 'first',
        'delta_mass': 'first',
        'spec_spec': 'first'
    }).reset_index()

    # Merge the counts back in
    df = agg_df.merge(group_counts, on=['ref_1_inchi', 'ref_2_inchi'])
    
    # for spec_spec, if False, then set ref_2_id to None
    df['ref_2_id'] = df.apply(lambda x: x['ref_2_id'] if x['spec_spec'] else None, axis=1)
    
    # drop cols
    df = df.drop(columns=['spec_spec']).reset_index(drop=True)
    
    # some more to deal with qry_id and qry_mz
    df['dataset'] = df['qry_id'].apply(lambda x: x.split(':')[0])
    df['file_path'] = df['qry_id'].apply(lambda x: x.split(':')[1])
    df['file_scan'] = df['qry_id'].apply(lambda x: x.split(':')[-1])
    
    # save as pkl
    df.to_pickle(out_path.replace('.parquet', '.pkl'))
    
    # save as parquet file
    df.to_parquet(out_path, index=False)  # qry_id, qry_mz, ref_1_id, ref_2_id, delta_mass, dataset, file_path, file_scan
    

def correct_qry_file_path(df_path, mri_path='/home/shipei/projects/revcos/search/uniquemri_refined.pkl', pos=True):
    
    df = pd.read_pickle(df_path)    
    df['short_mri'] = df['qry_id'].apply(lambda x: x.split(':scan')[0])
    
    print("Correcting qry file paths...")
    # load mri
    mri = pd.read_pickle(mri_path)
    
    # Create a mapping of ambiguous short_mri values
    print("Identifying ambiguous short_mri values...")
    ambiguous_short_mri = mri['short_mri'].value_counts()
    ambiguous_short_mri = set(ambiguous_short_mri[ambiguous_short_mri > 1].index)
    print(f"MRI dataframe: {len(mri)} rows")
    print(f"MRI dataframe: {mri['short_mri'].nunique()} unique short_mri values")
    print(f"MRI dataframe: {len(ambiguous_short_mri)} ambiguous short_mri values")
    
    # Split the dataframe into ambiguous and non-ambiguous rows
    non_ambiguous_df = df[~df['short_mri'].isin(ambiguous_short_mri)].copy()
    ambiguous_df = df[df['short_mri'].isin(ambiguous_short_mri)].copy()
    
    print(f"Processing {len(non_ambiguous_df)} non-ambiguous rows...")
    # For non-ambiguous rows, do a simple merge
    # Create a filtered mri dataframe with only the non-ambiguous short_mri values
    non_ambiguous_mri = mri.drop_duplicates(subset=['short_mri'])
    non_ambiguous_df = non_ambiguous_df.merge(
        non_ambiguous_mri[['short_mri', 'filepath']], 
        on='short_mri', 
        how='left'
    )
    non_ambiguous_df.rename(columns={'filepath': 'full_file_path'}, inplace=True)
    
    ambiguous_mri = mri[mri['short_mri'].isin(ambiguous_short_mri)].copy()
    print(f"{len(ambiguous_mri)} ambiguous rows in MRI dataframe")
    print(f"{ambiguous_mri['short_mri'].nunique()} unique short_mri values in ambiguous MRI dataframe")
    
    # Process ambiguous rows individually if there are any
    if len(ambiguous_df) > 0:
        print(f"Processing {len(ambiguous_df)} ambiguous rows...")
        ambiguous_df['full_file_path'] = None
        
        # Process each row individually instead of grouping by short_mri
        for idx, row in tqdm(ambiguous_df.iterrows(), total=len(ambiguous_df), desc="Matching ambiguous file paths"):
            short_mri = row['short_mri']
            scan_no = row['file_scan']
            qry_mz = row['qry_mz']
            
            # Get all matching rows from the MRI table
            if pos:
                matched_rows = ambiguous_mri[(ambiguous_mri['short_mri'] == short_mri) &
                                             (ambiguous_mri['polarity'] != 'neg')]
            else:
                matched_rows = ambiguous_mri[(ambiguous_mri['short_mri'] == short_mri) &
                                             (ambiguous_mri['polarity'] != 'pos')]
            
            if len(matched_rows) == 1:
                # If there's only one match, use it directly
                ambiguous_df.at[idx, 'full_file_path'] = matched_rows.iloc[0]['filepath']
                continue
            elif len(matched_rows) == 0:
                matched_rows = ambiguous_mri[(ambiguous_mri['short_mri'] == short_mri)]
            else:
                pass
            # Try to find a matching filepath
            for _, matched_row in matched_rows.iterrows():
                this_usi = matched_row['usi'] + ':scan:' + str(scan_no)
                prec_mz = load_from_usi(this_usi)
                    
                if prec_mz is not None and abs(prec_mz - qry_mz) < 1e-2:
                    # Update this specific row
                    ambiguous_df.at[idx, 'full_file_path'] = matched_row['filepath']
                    print(f"Matched {this_usi} with {matched_row['filepath']}")
                    break
            
            # If no match was found
            if ambiguous_df.at[idx, 'full_file_path'] is None:
                print(f"No match found for {short_mri} with scan {scan_no} and mz {qry_mz}")
        
        # Combine the results
        df = pd.concat([non_ambiguous_df, ambiguous_df], axis=0)
    else:
        df = non_ambiguous_df
    
    # Check for missing file paths
    missing_count = df['full_file_path'].isnull().sum()
    if missing_count > 0:
        print(f"Warning: {missing_count} file paths are still missing after matching.")
        missing_examples = df[df['full_file_path'].isnull()]['dataset'].tolist()
        missing_examples = set(missing_examples)
        print(f"Missing datasets: {missing_examples}")
    
    # drop cols
    df = df.drop(columns=['qry_id', 'qry_mz', 'file_path'])

    # save as pkl
    df.to_pickle(df_path.replace('.pkl', '_full_path.pkl'))
    
    # save as parquet file
    df.to_parquet(df_path.replace('.pkl', '_full_path.parquet'), index=False)
    return df

    
def add_search_results(result_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_best_raw.tsv', mode='pos',
                       out_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_web_2_full_path.parquet'):

    print(f"ADD search results from {result_path}")
    alldf = pd.read_csv(result_path, sep='\t', low_memory=False)
    
    print('Total rows:', len(alldf))
    
    alldf['ref_1_prec_int'] = alldf.apply(lambda x: max(x['ref_1_prec_frag_int'], x['ref_1_prec_frag_water_loss_int']), axis=1)
    
    # drop cols ref_1_prec_frag_int, ref_1_prec_frag_water_loss_int, ref_1_prec_nl_int, ref_1_prec_nl_water_loss_int, total_usage
    alldf = alldf.drop(columns=['ref_1_prec_frag_int', 'ref_1_prec_frag_water_loss_int', 'ref_1_prec_nl_int',
                               'ref_1_prec_nl_water_loss_int', 'total_usage']).reset_index(drop=True)
    
    alldf = alldf[(alldf['ref_1_score'] >= 0.8) & (alldf['ref_1_peak'] >= 3) & 
                  (alldf['ref_1_usage'] >= 0.30) & (alldf['ref_1_prec_int'] >= 0.05)].reset_index(drop=True)
    
    df_1 = alldf[(alldf['ref_2_score'] >= 0.6) & (alldf['ref_2_peak'] >= 3)].reset_index(drop=True)
    df_1['spec_spec'] = True
    
    df_2 = alldf[(alldf['ref_2_id'].isna())].reset_index(drop=True)
    df_2['spec_spec'] = False
    
    df = pd.concat([df_1, df_2], axis=0).reset_index(drop=True)
    del alldf
    print(f"Total rows after filtering: {len(df)}")
    
    # load metadata
    print("Loading metadata...")
    metadata = pd.read_pickle('/home/shipei/projects/revcos/search/all_ms2db_metadata.pkl')
    db_id_to_db = metadata.set_index('db_id')['db'].to_dict()
    db_id_to_mass = metadata.set_index('db_id')['monoisotopic_mass'].to_dict()
    db_id_to_inchi14 = metadata.set_index('db_id')['inchikey_14'].to_dict()
    
    # add metadata to df
    df['ref_1_mass'] = df['ref_1_id'].map(db_id_to_mass)
    df['delta_mass'] = df['qry_mz'] - df['ref_1_mass'] - 1.007276466812
    df['ref_1_db'] = df['ref_1_id'].map(db_id_to_db)
    df['ref_1_db'] = df['ref_1_db'].apply(lambda x: 1 if x == 'gnps' or x == 'mona' else 0)
    
    def _fill_ref_2_inchi(row):
        if pd.isnull(row['ref_2_id']):
            return str(round(row['delta_mass'], 2))
        else:
            id = str(row['ref_2_id'])
            if id.endswith('.0'):
                id = id[:-2]
            return db_id_to_inchi14[id]
    def _fill_ref_2_db(row):
        if pd.isnull(row['ref_2_id']):
            return 0
        else:
            id = str(row['ref_2_id'])
            if id.endswith('.0'):
                id = id[:-2]
            db = db_id_to_db[id]
            if db == 'gnps' or db == 'mona':
                return 1
            else:
                return 0
    df['ref_2_inchi'] = df.apply(_fill_ref_2_inchi, axis=1)
    df['ref_2_db'] = df.apply(_fill_ref_2_db, axis=1)
    
    df = df[df['delta_mass'] > 1.5].reset_index(drop=True)
    
    # sort by ref_1_inchi, ref_2_inchi, ref_1_db, ref_2_db, ref_1_prec_int, ref_1_usage, ref_1_score, ref_2_score
    df = df.sort_values(by=['ref_1_inchi', 'ref_2_inchi', 'ref_1_db', 'ref_2_db', 'ref_1_prec_int', 'ref_1_usage',
                            'ref_1_score', 'ref_2_score'], ascending=[True, True, False, False, False, False, False, False])
    df = df.reset_index(drop=True)
    
    # First get the group counts
    print("Grouping by ref_1_inchi and ref_2_inchi...")
    group_counts = df.groupby(['ref_1_inchi', 'ref_2_inchi']).size().reset_index(name='count')

    # Aggregation
    agg_df = df.groupby(['ref_1_inchi', 'ref_2_inchi']).agg({
        'qry_id': 'first',
        'qry_mz': 'first',
        'ref_1_id': 'first',
        'ref_2_id': 'first',
        'delta_mass': 'first',
        'spec_spec': 'first'
    }).reset_index()

    # Merge the counts back in
    df = agg_df.merge(group_counts, on=['ref_1_inchi', 'ref_2_inchi'])
    
    # for spec_spec, if False, then set ref_2_id to None
    df['ref_2_id'] = df.apply(lambda x: x['ref_2_id'] if x['spec_spec'] else None, axis=1)
    
    # drop cols
    df = df.drop(columns=['spec_spec']).reset_index(drop=True)
    
    # some more to deal with qry_id and qry_mz
    df['dataset'] = df['qry_id'].apply(lambda x: x.split(':')[0])
    df['file_path'] = df['qry_id'].apply(lambda x: x.split(':')[1])
    df['file_scan'] = df['qry_id'].apply(lambda x: x.split(':')[-1])
    
    #####################
    # correct qry file paths
    df['short_mri'] = df['qry_id'].apply(lambda x: x.split(':scan')[0])

    print("Loading existing qry file with full path...")        
    # first load the existing results
    calced_df = pd.read_parquet(f'/home/shipei/projects/revcos/search/results/{mode}_stitched/{mode}_web_2_full_path.parquet')
    
    def _get_merged_col(row):
        if pd.isna(row['ref_2_inchi']):
            return row['short_mri'] + '_' + row['ref_1_inchi'] + '_' + 'None'
        else:
            return row['short_mri'] + '_' + row['ref_1_inchi'] + '_' + row['ref_2_inchi']
    df['merged_col'] = df.apply(lambda x: _get_merged_col(x), axis=1)
    calced_df['merged_col'] = calced_df.apply(lambda x: _get_merged_col(x), axis=1)
    
    df = df.merge(calced_df[['merged_col', 'full_file_path']], on='merged_col', how='left')
    # drop merged_col
    df = df.drop(columns=['merged_col']).reset_index(drop=True)
    
    df_filled = df[df['full_file_path'].notna()].reset_index(drop=True).copy()
    df = df[df['full_file_path'].isna()].drop(columns=['full_file_path']).reset_index(drop=True)  # missing file paths
    
    # df_filled.to_parquet(out_path.replace('.parquet', '_filled.parquet'), index=False)
    #######
    # load mri
    mri = pd.read_pickle('/home/shipei/projects/revcos/search/uniquemri_refined.pkl')
    
    # Create a mapping of ambiguous short_mri values
    print("Identifying ambiguous short_mri values...")
    ambiguous_short_mri = mri['short_mri'].value_counts()
    ambiguous_short_mri = set(ambiguous_short_mri[ambiguous_short_mri > 1].index)
    print(f"MRI dataframe: {len(mri)} rows")
    print(f"MRI dataframe: {mri['short_mri'].nunique()} unique short_mri values")
    print(f"MRI dataframe: {len(ambiguous_short_mri)} ambiguous short_mri values")
    
    # Split the dataframe into ambiguous and non-ambiguous rows
    non_ambiguous_df = df[~df['short_mri'].isin(ambiguous_short_mri)].copy()
    ambiguous_df = df[df['short_mri'].isin(ambiguous_short_mri)].copy()
    
    print(f"Processing {len(non_ambiguous_df)} non-ambiguous rows...")
    # For non-ambiguous rows, do a simple merge
    # Create a filtered mri dataframe with only the non-ambiguous short_mri values
    non_ambiguous_mri = mri.drop_duplicates(subset=['short_mri'])
    non_ambiguous_df = non_ambiguous_df.merge(
        non_ambiguous_mri[['short_mri', 'filepath']], 
        on='short_mri', 
        how='left'
    )
    non_ambiguous_df.rename(columns={'filepath': 'full_file_path'}, inplace=True)
    
    ambiguous_mri = mri[mri['short_mri'].isin(ambiguous_short_mri)].copy()
    print(f"{len(ambiguous_mri)} ambiguous rows in MRI dataframe")
    print(f"{ambiguous_mri['short_mri'].nunique()} unique short_mri values in ambiguous MRI dataframe")
    
    # Process ambiguous rows individually if there are any
    if len(ambiguous_df) > 0:
        print(f"Processing {len(ambiguous_df)} ambiguous rows...")
        ambiguous_df['full_file_path'] = None
        
        # Process each row individually instead of grouping by short_mri
        for idx, row in tqdm(ambiguous_df.iterrows(), total=len(ambiguous_df), desc="Matching ambiguous file paths"):
            short_mri = row['short_mri']
            scan_no = row['file_scan']
            qry_mz = row['qry_mz']
            
            # Get all matching rows from the MRI table
            if mode == 'pos':
                matched_rows = ambiguous_mri[(ambiguous_mri['short_mri'] == short_mri) &
                                             (ambiguous_mri['polarity'] != 'neg')]
            else:
                matched_rows = ambiguous_mri[(ambiguous_mri['short_mri'] == short_mri) &
                                             (ambiguous_mri['polarity'] != 'pos')]
            
            if len(matched_rows) == 1:
                # If there's only one match, use it directly
                ambiguous_df.at[idx, 'full_file_path'] = matched_rows.iloc[0]['filepath']
                continue
            elif len(matched_rows) == 0:
                matched_rows = ambiguous_mri[(ambiguous_mri['short_mri'] == short_mri)]
            else:
                pass
            # Try to find a matching filepath
            for _, matched_row in matched_rows.iterrows():
                this_usi = matched_row['usi'] + ':scan:' + str(scan_no)
                prec_mz = load_from_usi(this_usi)
                    
                if prec_mz is not None and abs(prec_mz - qry_mz) < 1e-2:
                    # Update this specific row
                    ambiguous_df.at[idx, 'full_file_path'] = matched_row['filepath']
                    # print(f"Matched {this_usi} with {matched_row['filepath']}")
                    break
            
            # If no match was found
            if ambiguous_df.at[idx, 'full_file_path'] is None:
                print(f"No match found for {short_mri} with scan {scan_no} and mz {qry_mz}")
        
        # Combine the results - FIX: Reset index on both dataframes before concatenation
        # non_ambiguous_df.to_parquet(out_path.replace('.parquet', '_non_ambiguous.parquet'), index=False)
        # ambiguous_df.to_parquet(out_path.replace('.parquet', '_ambiguous.parquet'), index=False)
        df = pd.concat([non_ambiguous_df, ambiguous_df], axis=0, ignore_index=True).reset_index(drop=True)
    else:
        df = non_ambiguous_df
    
    # Check for missing file paths
    missing_count = df['full_file_path'].isnull().sum()
    if missing_count > 0:
        print(f"Warning: {missing_count} file paths are still missing after matching.")
        missing_examples = df[df['full_file_path'].isnull()]['dataset'].tolist()
        missing_examples = set(missing_examples)
        print(f"Missing datasets: {missing_examples}")
    
    # concat
    # df.to_parquet(out_path.replace('.parquet', '_missing.parquet'), index=False)
    df = pd.concat([df_filled, df], axis=0, ignore_index=True).reset_index(drop=True)
    # drop cols
    df = df.drop(columns=['qry_id', 'qry_mz', 'file_path'])

    # save as pkl
    df.to_pickle(out_path.replace('.parquet', '.pkl'))
    
    # save as parquet file
    df.to_parquet(out_path, index=False)
    return df
    


from requests import get
from json import loads
from time import sleep


def load_from_usi(usi, max_retries=3, delay=1):
    for attempt in range(max_retries):
        try:
            url = 'https://metabolomics-usi.gnps2.org/json/?usi1=' + usi
            response = get(url, timeout=5)
            json_data = loads(response.text)

            # check if the USI is valid
            if 'error' in json_data:
                return None

            prec_mz = json_data['precursor_mz']

            return float(prec_mz)
        except:
            if attempt < max_retries - 1:  # Don't sleep on the last attempt
                sleep(delay)
            continue

    return None


def prepare_ms2db_metadata(path='/Users/shipei/Documents/projects/conjugated_metabolome/db/ms2db/all/all_ms2db_metadata.pkl'):
    
    metadata = pd.read_pickle(path)
    
    metadata = metadata[['name', 'db', 'db_id', 'inchikey_14', 'monoisotopic_mass']]
    # change the column 'db' to categorical
    metadata['db'] = metadata['db'].astype('category')
    # change the column 'monoisotopic_mass' to 4 decimal places
    metadata['monoisotopic_mass'] = metadata['monoisotopic_mass'].round(4)
    
    # show the first 5 rows of gnps db
    print(metadata['db'].cat.categories)
    # for each category, show the first 5 rows
    for category in metadata['db'].cat.categories:
        print(f"Category: {category}")
        print(metadata[metadata['db'] == category].head(5))

    # verify if all gnps spectra, db_id starts with 'CCMSLIB000'
    print(metadata[metadata['db'] == 'gnps']['db_id'].str.startswith('CCMSLIB000').all())
    # verify if all massbank spectra, db_id starts with 'MSBNK-'
    print(metadata[metadata['db'] == 'massbank']['db_id'].str.startswith('MSBNK-').all())
    
    # some rules
    # 1. if db is 'gnps', then db_id removes 'CCMSLIB000'
    # 2. if db is 'massbank', then db_id removes 'MSBNK-'
    # 3. if db is mona or nist20, then db_id remains the same
    # def remove_prefix(row):
    #     if row['db'] == 'gnps':
    #         return row['db_id'].replace('CCMSLIB000', '')
    #     elif row['db'] == 'massbank':
    #         return row['db_id'].replace('MSBNK-', '')
    #     else:
    #         return row['db_id']
    def remove_prefix(db_id):
        if pd.isnull(db_id):
            return db_id
        return db_id.replace('CCMSLIB000', '').replace('MSBNK-', '')
    metadata['db_id'] = metadata['db_id'].apply(remove_prefix)
    
    # check if all db_id are unique
    if metadata['db_id'].is_unique:
        print("All db_id are unique")
    
    # # drop col 'db'
    # metadata = metadata.drop(columns=['db'])
    
    # save as parquet file
    metadata.to_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/ms2db.parquet', index=False)


def prepare_mri_table():
    mri = pd.read_csv('/Users/shipei/Documents/projects/conjugated_metabolome/db/mri/uniquemri.csv', low_memory=False)
    mri = mri[['dataset', 'filepath', 'usi', 'spectra_ms1', 'spectra_ms2']]
    
    mri['file_type'] = mri['filepath'].apply(lambda x: os.path.splitext(x)[1].replace('.', '').lower())
    print(mri['file_type'].value_counts())
    mri = mri[mri['file_type'].isin(['mzml', 'mzxml'])].reset_index(drop=True)
    
    matched_datasets = pd.read_pickle('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/matched_datasets.pkl')
    mri = mri[mri['dataset'].isin(matched_datasets)].reset_index(drop=True)
    print(f"Total number of rows after filtering: {len(mri)}")
    
    # Create short_mri and polarity columns
    mri['short_mri'] = mri.apply(lambda x: x['dataset'] + ':' + os.path.basename(os.path.splitext(x['filepath'])[0]), axis=1)
    mri['polarity'] = mri['filepath'].apply(lambda x: determine_pos_neg(x))
    
    mri_1 = mri[(mri['spectra_ms1'] > 0) | (mri['spectra_ms2'] > 0)].reset_index(drop=True)
    # deduplicate mri_1
    mri_1 = mri_1.drop_duplicates(subset=['short_mri', 'polarity', 'spectra_ms1', 'spectra_ms2'], keep='first').reset_index(drop=True)
    mri_2 = mri[(mri['spectra_ms1'] == 0) & (mri['spectra_ms2'] == 0)].reset_index(drop=True)
    mri = pd.concat([mri_1, mri_2], axis=0).reset_index(drop=True)
    
    # For short_mri + polarity groups with exactly two rows (mzml and mzxml),
    # if one has spectra_ms2 > 0 and the other has spectra_ms2 = 0, keep the one with spectra_ms2 > 0
    print("Identifying groups with exactly one mzml and one mzxml file...")

    # First, identify all the short_mri + polarity groups that have exactly 2 rows
    group_sizes = mri.groupby(['short_mri', 'polarity']).size()
    exact_two_groups = group_sizes[group_sizes == 2].index.tolist()
    print(f"Found {len(exact_two_groups)} short_mri + polarity groups with exactly 2 rows")

    # Filter for just those groups, much faster than processing everything
    exact_two_df = mri[mri.set_index(['short_mri', 'polarity']).index.isin(exact_two_groups)]

    # Use vectorized operations to identify the groups with one mzml and one mzxml
    format_counts = exact_two_df.groupby(['short_mri', 'polarity'])['file_type'].apply(lambda x: set(x))
    mixed_format_groups = format_counts[format_counts == {'mzml', 'mzxml'}].index.tolist()
    print(f"Found {len(mixed_format_groups)} short_mri + polarity groups with one mzml and one mzxml file")

    if len(mixed_format_groups) > 0:
        # Filter down to just the mixed format groups
        mixed_format_df = exact_two_df[exact_two_df.set_index(['short_mri', 'polarity']).index.isin(mixed_format_groups)]
        
        # Create a mask for rows to remove - where spectra_ms2 = 0 and there's an alternative with spectra_ms2 > 0
        # This is done in one step without looping
        mixed_format_groups_with_inconsistent_ms2 = mri.loc[mixed_format_df.index].groupby(['short_mri', 'polarity'])['spectra_ms2'].agg(['max', 'min'])
        inconsistent_groups = mixed_format_groups_with_inconsistent_ms2[(mixed_format_groups_with_inconsistent_ms2['max'] > 0) & 
                                                                    (mixed_format_groups_with_inconsistent_ms2['min'] == 0)].index
        
        if len(inconsistent_groups) > 0:
            print(f"Found {len(inconsistent_groups)} groups where one file has MS2 spectra and the other doesn't")
            # Identify rows to remove
            rows_to_remove = mixed_format_df[
                mixed_format_df.set_index(['short_mri', 'polarity']).index.isin(inconsistent_groups) & 
                (mixed_format_df['spectra_ms2'] == 0)
            ].index
            
            print(f"Removing {len(rows_to_remove)} rows with zero MS2 spectra where alternatives exist")
            mri = mri.drop(index=rows_to_remove).reset_index(drop=True)
        else:
            print("No inconsistent MS2 spectra found in mixed format groups")
    else:
        print("No mixed format groups found")
    
    print("Prioritizing files based on spectra counts...")
    # Sort the dataframe to prioritize files with more spectra
    mri = mri.sort_values(
        by=['short_mri', 'spectra_ms1', 'spectra_ms2', 'file_type'],
        ascending=[True, False, False, True]
    )
    
    # ccms_peak and raw files are not needed
    mri['std_usi'] = mri['usi'].apply(lambda x: x.replace('ccms_peak/', 'peak/').replace('peaks/', 'peak/').replace('raw/', 'peak/').replace('mzml/', 'peak/').replace('mzML/', 'peak/').replace('mzXML/', 'peak/').replace('mzxml/', 'peak/').replace('mzmls/', 'peak/').replace('mzMLs/', 'peak/').replace('mzXMLs/', 'peak/').replace('mzxmls/', 'peak/'))
    mri = mri.drop_duplicates(subset=['std_usi'], keep='first').reset_index(drop=True)
    mri = mri.drop(columns=['std_usi'])
    
    print(f"Total number of rows after deduplication: {len(mri)}")
    print(f"Total number of unique short_mri: {mri['short_mri'].nunique()}")
   
    # print the number of rows with the same short_mri and polarity
    print("Checking for duplicate short_mri + polarity combinations...")
    duplicates = mri.groupby(['short_mri', 'polarity']).size().reset_index(name='count')
    duplicates = duplicates[duplicates['count'] > 1]
    if len(duplicates) > 0:
        print(f"Found {len(duplicates)} duplicate short_mri + polarity combinations")
        print(duplicates)
    else:
        print("No duplicate short_mri + polarity combinations found")
    
    # # Then keep only the first occurrence of each short_mri + polarity combination
    # mri = mri.drop_duplicates(subset=['short_mri', 'polarity', 'spectra_ms1', 'spectra_ms2'], keep='first').reset_index(drop=True)
    
    print(f"Total number of rows: {len(mri)}")
    print(f"Total number of unique short_mri: {mri['short_mri'].nunique()}")
    print(f"Total number of unique short_mri + polarity combinations: {mri.groupby(['short_mri', 'polarity']).ngroups}")
    
    # save as tsv and pkl
    mri.to_csv('/Users/shipei/Documents/projects/conjugated_metabolome/db/mri/uniquemri_refined.tsv', sep='\t', index=False)
    mri.to_pickle('/Users/shipei/Documents/projects/conjugated_metabolome/db/mri/uniquemri_refined.pkl')


def determine_pos_neg(text):   
   lower_text = str(text).lower()
   pos = lower_text.count('pos')
   neg = lower_text.count('neg')
   
   if pos > neg:
       return 'pos'
   elif neg > pos:
       return 'neg'
   else:
       return 'unknown'


def prepare_matched_mri_list():
    neg = pd.read_csv('/home/shipei/projects/revcos/search/results/neg_stitched/neg_best_raw.tsv', sep='\t')
    neg['dataset'] = neg['qry_id'].apply(lambda x: x.split(':')[0])
    neg_datasets = neg['dataset'].unique().tolist()
    del neg
    
    pos = pd.read_csv('/home/shipei/projects/revcos/search/results/pos_stitched/pos_best_raw.tsv', sep='\t')
    pos['dataset'] = pos['qry_id'].apply(lambda x: x.split(':')[0])
    pos_datasets = pos['dataset'].unique().tolist()
    del pos
    
    all_datasets = list(set(neg_datasets + pos_datasets))
    
    # save as pkl
    with open('/home/shipei/projects/revcos/search/matched_datasets.pkl', 'wb') as f:
        pickle.dump(all_datasets, f)
    print(f"Saved {len(all_datasets)} datasets to matched_datasets.pkl")
    

def refine_search_results(result_path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/neg_web_full_path.parquet',
                          out_path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/neg_refined.parquet'):
        
    print(f"Refining search results from {result_path} to {out_path}")
    df = pd.read_parquet(result_path)
    
    # load mri
    mri = pd.read_csv('/Users/shipei/Documents/projects/conjugated_metabolome/db/mri/uniquemri.csv', low_memory=False)
    mri = mri[['dataset', 'filepath', 'usi', 'spectra_ms1', 'spectra_ms2']]        
    mri['file_type'] = mri['filepath'].apply(lambda x: os.path.splitext(x)[1].replace('.', '').lower())
    mri = mri[mri['file_type'].isin(['mzml', 'mzxml'])].reset_index(drop=True)
    
    matched_datasets = pd.read_pickle('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/matched_datasets.pkl')
    mri = mri[mri['dataset'].isin(matched_datasets)].reset_index(drop=True)
    mri['short_mri'] = mri.apply(lambda x: x['dataset'] + ':' + os.path.basename(os.path.splitext(x['filepath'])[0]), axis=1)
    
    print('None full_file_path:', df['full_file_path'].isnull().sum())
    
    df_missing = df[df['full_file_path'].isnull()].reset_index(drop=True).copy()
    df_filled = df[df['full_file_path'].notnull()].reset_index(drop=True).copy()
    
    for i, row in tqdm(df_missing.iterrows(), total=len(df_missing), desc="Processing rows"):
        if not pd.isnull(row['full_file_path']):
            continue
        # get matches by short_mri
        short_mri = row['short_mri']
        matches = mri[mri['short_mri'] == short_mri]
        if matches.empty:
            continue
        # if there are multiple matches, choose the first one
        df_missing.at[i, 'full_file_path'] = matches.iloc[0]['filepath']
        
    print('None full_file_path after processing:', df['full_file_path'].isnull().sum())
    
    # merge df_filled and df_missing
    df = pd.concat([df_filled, df_missing], axis=0, ignore_index=True).reset_index(drop=True)
    
    # delta_mass, round to 2 decimal places
    df['delta_mass'] = df['delta_mass'].round(2)
    
    def remove_prefix(db_id):
        if pd.isnull(db_id):
            return db_id
        return db_id.replace('CCMSLIB000', '').replace('MSBNK-', '')
    df['ref_1_id'] = df['ref_1_id'].apply(remove_prefix)
    df['ref_2_id'] = df['ref_2_id'].apply(remove_prefix)
    
    # drop cols of ref_1_inchi, ref_2_inchi, short_mri
    df = df.drop(columns=['ref_1_inchi', 'ref_2_inchi', 'short_mri']).reset_index(drop=True)
    # # save as pkl
    # df.to_pickle(out_path.replace('.parquet', '.pkl'))
    
    # save as parquet file
    df.to_parquet(out_path, index=False)


def compress_ms2db_metadata(path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/ms2db.parquet'):
    """
    Compress the ms2db metadata file, by db ids that are in the result files, to reduce its size.
    """
    db = pd.read_parquet(path)

    # load the search results
    print("Loading search results to filter db_ids...")
    neg = pd.read_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/neg_refined.parquet')
    pos = pd.read_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/pos_refined.parquet')
    all_ids = set(neg['ref_1_id']).union(set(neg['ref_2_id'])).union(set(pos['ref_1_id'])).union(set(pos['ref_2_id']))
    del neg, pos
    
    print(f"Total unique db_ids in search results: {len(all_ids)}")
    # filter the db by db_id
    db = db[db['db_id'].isin(all_ids)].reset_index(drop=True)
    print(f"Total unique db_ids in ms2db metadata after filtering: {len(db)}")
    # save as parquet file
    db.to_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/ms2db.parquet', index=False)


def refine_by_count():
    df = pd.read_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/neg_refined.parquet')
    df = df[df['count'] >= 3].reset_index(drop=True)
    df.to_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/neg_refined.parquet', index=False)
    
    df = pd.read_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/pos_refined.parquet')
    df = df[df['count'] >= 3].reset_index(drop=True)
    df.to_parquet('/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/pos_refined.parquet', index=False)


if __name__ == "__main__":
    import os
    os.chdir(os.path.dirname(__file__))
     
    # prepare_ms2db_metadata() # local
    
    # prepare_matched_mri_list()  # server
    # prepare_mri_table()  # local
    
    ######first batch, prec int >= 0.1
    # server
    # prepare_search_results(result_path='/home/shipei/projects/revcos/search/results/neg_stitched/neg_best_raw.tsv',
    #                        out_path='/home/shipei/projects/revcos/search/results/neg_stitched/neg_web.parquet')
    
    # prepare_search_results(result_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_best_raw.tsv',
    #                        out_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_web.parquet')

    # correct_qry_file_path('/home/shipei/projects/revcos/search/results/neg_stitched/neg_web.pkl',
    #                       mri_path='/home/shipei/projects/revcos/search/uniquemri_refined.pkl', pos=False)
    # correct_qry_file_path('/home/shipei/projects/revcos/search/results/pos_stitched/pos_web.pkl',
    #                       mri_path='/home/shipei/projects/revcos/search/uniquemri_refined.pkl', pos=True)
        
    ######2nd batch, prec int >= 0.05
    # add_search_results(result_path='/home/shipei/projects/revcos/search/results/neg_stitched/neg_best_raw.tsv', mode='neg',
    #                    out_path='/home/shipei/projects/revcos/search/results/neg_stitched/neg_web_2_full_path.parquet')
    # add_search_results(result_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_best_raw.tsv', mode='pos',
    #                    out_path='/home/shipei/projects/revcos/search/results/pos_stitched/pos_web_2_full_path.parquet')
    
    ################
    # refine_search_results(result_path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/neg_web_2_full_path.parquet',
    #                       out_path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/neg_refined.parquet')
    
    # refine_search_results(result_path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/pos_web_2_full_path.parquet',
    #                       out_path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/pos_refined.parquet')
    
    ##############
    # refine by count
    refine_by_count()
    
    ##############    
    compress_ms2db_metadata(path='/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/ms2db.parquet')