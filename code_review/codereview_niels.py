import pandas as pd
import numpy as np

def go_terms(id_list, go_dict):
    """
    Takes the UniProt ID column (as a list) and the dictionary translating
    IDs to corresponding GO terms. Returns a list that can be used as a new
    column in the dataframe containing all GO-terms annotated to that protein
    as a set.
    """
    return [go_dict[x] if x in go_dict.keys() else np.nan for x in id_list]

def interesting_gos(go_interesting_file):
    """
    Takes a file containing a GO-term on each line and returns a set
    of those terms.
    """
    with open(go_interesting_file) as file:
        return set([x.strip() for x in file])

def go_terms_hits(all_terms_column, go_interesting):
    """
    Takes a column containing sets of GO-terms and returns a
    list that can be used as a new column containing the GO-terms that are
    in the interesting GO terms set.
    """
    all_terms_column = list(all_terms_column)
    interesting_annotations_column = []
    for annotations in all_terms_column:
        newset = set()
        if isinstance(annotations, set):
            for annotation in annotations:
                if annotation in go_interesting:
                    newset.add(annotation)
        if len(newset) == 0:
            interesting_annotations_column.append(np.nan)
        else:
            interesting_annotations_column.append(newset)
    return interesting_annotations_column

def interesting_go_counts(df, go_interesting, branch):
    """
    Takes the complete cofactorslist, the set of interesting GO terms
    and a branch of the gene ontology ('MF', 'CC' or 'BP').
    Returns the dataframe with two new columns. One with the total number
    of terms annotated to each protein, one with the number of interesting
    annotations for each protein.
    """
    hits_col = []
    tot_ann_col = []
    fraction_col = []
    for index, row in list_df.iterrows():
        gos = row['{} terms'.format(branch)]
        hits = 0
        if isinstance(gos, set):
            annots = len(gos)
            for term in gos:
                if term in go_interesting:
                    hits += 1
        else:
            annots = 0
        hits_col.append(hits)
        tot_ann_col.append(annots)
    outdf = pd.DataFrame({
        '{} hits'.format(branch) : hits_col,
        '{} total annotations'.format(branch) : tot_ann_col,
    })
    return pd.concat([df, outdf], axis = 1)

def domain_scores_col(df, weightdict, branch):
    """
    Takes the complete cofactorlist, the dictionary with the weight
    that was assigned to each interesting GO-term and a branch of the
    GO-ontology.
    Returns the cofactorlist with two new columns. One column with the
    summed weight of the 5 GO terms with the highest weight for each protein,
    one column with a penalty score for proteins that do not have any
    interesting terms, but do have other terms.
    """
    new_col_pos = []
    new_col_neg = []
    for index, row in df.iterrows():
        terms = row['{} interesting'.format(branch)]
        total = row['{} total annotations'.format(branch)]
        if isinstance(terms, set):
            neg_score = 0
            all_weights = [weightdict[x] for x in terms]
            if len(all_weights) <= 5:
                weightsum = sum(all_weights)
            else:
                all_weights.sort(reverse=True)
                weightsum = sum(all_weights[:5])
        else:
            weightsum = 0
            if total > 5:
                neg_score = -5
            else:
                neg_score = -total
        new_col_pos.append(weightsum)
        new_col_neg.append(neg_score)
    df['{} go score'.format(branch)] = new_col_pos
    df['{} go score penalty'.format(branch)] = new_col_neg
    return df

def mine_score(df):
    """
    Takes the cofactorlist, returns the cofactorlist with a column containing
    the mining scores for each protein. (One point for each dataset it was
    found in).
    """
    in_humap = df['ComplexID_huMAP'].notnull().astype('int')
    in_corum = df['ComplexID_CORUM'].notnull().astype('int')
    in_biogrid = df['Interaction_TF_BioGRID']
    in_intact = df['Interaction_TF_IntACT']
    in_bait = (df['In_bait_crems'] |
               df['In_bait_snfs'] |
               df['In_bait_nursa'] |
               df['In_bait_gocofs'] |
               df['In_NVS'])

    score = (in_humap +
             in_corum +
             in_biogrid +
             in_intact +
             in_bait)
    
    df['Mine score'] = score
    return df

# Read in the complete cofactorlist.
list_df = pd.read_csv(cof_list)


# Add columns with all GO-terms annotated
# to each protein for each aspect of the GO.
unips = list(list_df['UniProt ID'])
list_df['CC terms'] = go_terms(unips, cc_dict)
list_df['BP terms'] = go_terms(unips, bp_dict)
list_df['MF terms'] = go_terms(unips, mf_dict)

# Create sets of interesting GO-terms in each
# aspect of the GO.
interesting_cc, interesting_bp, interesting_mf = (interesting_gos(x) for x in [go_interesting_cc_file,
                                                                               go_interesting_bp_file,
                                                                               go_interesting_mf_file])

# Add columns with all interesting GO-terms
# annotated to each protein for each aspect
# of the GO.
list_df['CC interesting'] = go_terms_hits(list_df['CC terms'], interesting_cc)
list_df['BP interesting'] = go_terms_hits(list_df['BP terms'], interesting_bp)
list_df['MF interesting'] = go_terms_hits(list_df['MF terms'], interesting_mf)

# Add columns with the counts on number of interesting
# annotations and total number of annotations for each
# protein for each aspect of the GO.
list_df = interesting_go_counts(list_df, interesting_cc, 'CC')
list_df = interesting_go_counts(list_df, interesting_bp, 'BP')
list_df = interesting_go_counts(list_df, interesting_mf, 'MF')

# Read in the text file with the interesting GO-terms
# and their scores. Change that score into a weight for
# each term and create a dictionary with the terms as keys
# and their weight as values.
weighting_df = pd.read_table('/mbshome/nvelthuijs/Cofactors/List2/Datafiles/weighting_goterms.txt',
                            names = ['Term', 'Counts'])
def calc_weight(count):
    return (count/587)/2+ 0.5
weighting_df['Weight'] = [calc_weight(row['Counts']) for index, row in weighting_df.iterrows()]
weighted_terms_dict = {row['Term']:row['Weight'] for index, row in weighting_df.iterrows()}

# Add columns with the scores and penalty for each protein in
# each aspect of the GO.
for x in ['CC', 'BP', 'MF']:
    list_df = domain_scores_col(list_df, weighted_terms_dict, x)

# Sum the GO-scores and -penalties to the total GO-score.
list_df['Total GO score'] = list_df['CC go score'] + list_df['BP go score'] + list_df['MF go score']
list_df['Total GO score penalty'] = list_df['CC go score penalty'] + list_df['BP go score penalty'] + list_df['MF go score penalty']
list_df['GO score'] = (list_df['Total GO score'] + list_df['Total GO score penalty']) / 3

# Add a column with the mining score.
list_df = mine_score(list_df)

# Sum the GO-score and mining score.
list_df['Final score'] = list_df['GO score'] + list_df['Mine score']

# Write to csv.
list_df.to_csv('List_scores.csv', index = False)
