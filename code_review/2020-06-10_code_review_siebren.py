# The seq2science tool allows for differential gene expression analysis.
# To perform this, the user specifies one or more designs in the `config.yaml` like so:
#diffexp:
#  deseq2:
#    contrasts:
#      - 'disease_control_disease'
#      - 'stage_2_1'
#      - 'stage_all_1'
#
# these must be checked (the columns must be present in the `samples.tsv` file, and the factors must occur within the columns)
# the designs must be parsed (batch effects and specific modules to be extracted)
# finally, Snakemake must be instructed which of these to perform.
# In true snakemake fashion, this happens in inverted order.

# this code was written quite a while ago. So it could most likely use some cleaning up :)


## in the Snakefile we specify which files need to be produced

contrasts = get_contrasts()

rule count_all:
    input:
        [expand(f"{{dge_dir}}/{{assemblies}}-{contrast}.diffexp.tsv", **{**config, **{'assemblies': set(samples['assembly'])}}) for contrast in contrasts] if config.get('contrasts', False) else [],


## Extracting contrasts from the config dict

def parse_DE_contrasts(de_contrast):
    """
    Extract batch and contrast groups from a DE contrast design
    """
    original_contrast = de_contrast

    # remove whitespaces (and '~'s if used)
    de_contrast = de_contrast.replace(" ", "").replace("~", "")

    # split and store batch effect
    batch = None
    if '+' in de_contrast:
        batch =  de_contrast.split('+')[0]
        de_contrast = de_contrast.split('+')[1]

    # parse contrast
    parsed_contrast = de_contrast.split('_')
    return original_contrast, parsed_contrast, batch


def get_contrasts():
    """
    splits contrasts that contain multiple comparisons
    """
    if not config.get('contrasts', False):
        return []

    # contrasts from config
    old_contrasts = list(config["contrasts"])

    new_contrasts = []
    for contrast in old_contrasts:
        original_contrast, contrast, batch = parse_DE_contrasts(contrast)

        l = len(contrast)
        if l == 1:
            # write out the full contrast
            lvls = samples[contrast[0]].dropna().unique().tolist()
            new_contrast = batch + '+' + contrast[0] + '_' + lvls[0] + '_' + lvls[1] if batch is not None else contrast[0] + '_' + lvls[0] + '_' + lvls[1]
            new_contrasts.append(new_contrast)
        elif l == 2 or (l == 3 and 'all' in contrast[1:]):
            # create a list of all contrasts designed by the 'groupA vs all' design
            reflvl = str(contrast[2]) if contrast[1] == 'all' else str(contrast[1])
            lvls = samples[contrast[0]].dropna().unique().tolist()
            lvls = list(map(lambda x: str(x), lvls))
            lvls.remove(reflvl)

            for lvl in lvls:
                new_contrast = batch + '+' + contrast[0] + '_' + lvl + '_' + reflvl if batch is not None else contrast[0] + '_' + lvl + '_'  + reflvl
                new_contrasts.append(new_contrast)
        else:
            # remove '~', for uniformity
            new_contrast = original_contrast.replace("~", "").replace(" ", "")
            new_contrasts.append(new_contrast)

    # get unique elements
    new_contrasts = list(set(new_contrasts))
    return new_contrasts


## sanity checks in the configuration step

if config.get('contrasts', False):
    # check differential gene expression contrasts
    old_contrasts = list(config["contrasts"])
    for contrast in old_contrasts:
        original_contrast, parsed_contrast, batch = parse_DE_contrasts(contrast)

        # Check if the column names can be recognized in the contrast
        assert parsed_contrast[0] in samples.columns and parsed_contrast[0] not in ["sample", "assembly"], \
            (f'\nIn contrast design "{original_contrast}", "{parsed_contrast[0]} ' +
             f'does not match any valid column name in {config["samples"]}.\n')
        if batch is not None:
            assert batch in samples.columns and batch not in ["sample", "assembly"], \
                (f'\nIn contrast design "{original_contrast}", the batch effect "{batch}" ' +
                 f'does not match any valid column name in {config["samples"]}.\n')

        # Check if the groups described by the contrast can be identified and found in samples.tsv
        l = len(parsed_contrast)
        assert l < 4, ("\nA differential expression contrast couldn't be parsed correctly.\n" +
                       f"{str(l-1)} groups were found in '{original_contrast}' " +
                       f"(groups: {', '.join(parsed_contrast[1:])}).\n\n" +
                       f'Tip: do not use underscores in the columns of {config["samples"]} referenced by your contrast.\n')
        if l == 1:
            # check if contrast column has exactly 2 factor levels (per assembly)
            tmp = samples[['assembly', parsed_contrast[0]]].dropna()
            factors = pd.DataFrame(tmp.groupby('assembly')[parsed_contrast[0]].nunique())
            assert all(factors[parsed_contrast[0]] == 2),\
                ('\nYour contrast design, ' + original_contrast +
                 ', contains only a column name (' + parsed_contrast[0] +
                 '). \nIf you wish to compare all groups in this column, add a reference group. ' +
                 'Number of groups found (per assembly): \n' + str(factors[parsed_contrast[0]]))
        else:
            # check if contrast column contains the groups
            for group in parsed_contrast[1:]:
                if group != 'all':
                    assert str(group) in [str(i) for i in samples[parsed_contrast[0]].tolist()],\
                        ('\nYour contrast design contains group ' + group +
                        ' which cannot be found in column ' + parsed_contrast[0] +
                         ' of ' + config["samples"] + '.\n')
