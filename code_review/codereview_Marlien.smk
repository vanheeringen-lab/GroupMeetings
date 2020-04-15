configfile: "config.yml"

import pandas as pd
samples = pd.read_table(config['samples'], sep="\t")
assemblies = set(samples["assembly"])


config['top_N'] = str(config['top_N'])
config['slop'] = str(config['slop'])
config['clusters'] = str(config['clusters'])

rule all:
    input:
        [config['out_dir'] + assembly + "_maelstrom.out/" for assembly in assemblies],
        [config['out_dir'] + "heatmaps/" + assembly for assembly in assemblies]

# def get_inputfiles
def get_inputfiles(wildcards):
    gsms = samples[samples["assembly"] == wildcards.species]["sample"]
    paths = []
    for gsm in gsms:
        paths.append(config['in_dir'] + "summits/" + wildcards.species + "-" + gsm + "_summits.bed")
    return paths

def get_inputfiles_bam(wildcards):
    gsms = samples[samples["assembly"] == wildcards.species]["sample"]
    paths = []
    for gsm in gsms:
        paths.append(config['in_dir'] + "bam/" + wildcards.species + "-" + gsm + ".samtools-coordinate.bam")
    return paths

rule combine_peaks:
    input:
        summitfiles = get_inputfiles,
        genome = config['in_dir'] + "genomes/{species}.fa"
    output:
        config['out_dir'] + "combinepeaks/{species}.bed"
    shell:
        "combine_peaks -i {input.summitfiles} -g {input.genome} > {output}"

rule bedtools_slop:
    input:
        bedfile = rules.combine_peaks.output,
        chromsize = config['in_dir'] + "{species}.chrom.sizes"
    params:
        width = config['slop']
    output:
        config['out_dir'] + "combinepeaks/" + "{species}_" + config['slop'] + ".bed"
    shell:
        "bedtools slop -i {input.bedfile} -g {input.chromsize} -b {params.width} > {output}"

rule coverage_table:
    input:
        bedfile = config['in_dir'] + "combinepeaks/{species}_" + config['slop'] + ".bed",
        bamfiles = get_inputfiles_bam
    output:
        config['out_dir'] + "{species}_cp" + config['slop'] + "_ct.txt"
    shell:
        "coverage_table -p {input.bedfile} -d {input.bamfiles} > {output}"

rule change_header:
    input:
        rules.coverage_table.output,
    output:
        config['out_dir'] + "{species}_stage_ct.txt"
    run:
        df = pd.read_csv(str(input), comment="#", index_col=0, sep="\t")

        header = list(df.columns.values)
        sample_stage = dict(zip(samples['sample'], samples['stage']))
        print(header)
        print(sample_stage)

        for sample, stage in sample_stage.items():
            for name in header:
                if sample in name:
                    df.columns = df.columns.str.replace(sample, stage)
                    df.columns = df.columns.str.replace(".samtools-coordinate", "")

        print(df.columns.values)
        df.to_csv(str(output), index = True, index_label = str, header=True, sep="\t")

rule quantile_normalization:
    input:
        rules.change_header.output
    output:
        config['out_dir'] + "qNormalization/{species}_ct_qN.txt"
    run:
        import pandas as pd
        import numpy as np
        def quantileNormalize(df_input):

            df = df_input.copy()

            # ranking the values and sort median on axis=1
            dic = {}
            for col in df:
                dic.update({col : sorted(df[col])})
            sorted_df = pd.DataFrame(dic)
            rank = sorted_df.median(axis = 1).tolist()

            # sort
            for col in df:
                t = np.searchsorted(np.sort(df[col]), df[col])
                df[col] = [rank[i] for i in t]

            print(df)
            return df

        species = pd.read_csv(str(input), comment='#', index_col=0, sep="\t")
        species_qN = quantileNormalize(species)
        species_qN.to_csv(str(output), index = True, index_label = str, header=True, sep="\t")

rule transform_log_top_mc:
    input:
        rules.quantile_normalization.output
    params:
        toppeaks = int(config['top_N'])
    output:
        output_csv = config['out_dir'] + "{species}_log_diftop.txt",
        output_bed = config['out_dir'] + "{species}_log_diftop.bed"
    run:
        import pandas as pd
        import numpy as np
        from sklearn import preprocessing

        # import coveragte table
        input_ct = pd.read_csv(str(input), comment="#", index_col=0, sep="\t")

        # log transformation
        log_table = np.log1p(input_ct)

        # log_table max-min in a row
        df_dif = log_table.max(axis=1) - log_table.min(axis=1)
        top = df_dif.sort_values().tail(params.toppeaks).index
        output_top = log_table.loc[top]

        # centered on zero
        mc = preprocessing.scale(output_top, axis=1, with_mean='True', with_std=False)
        output_mc = pd.DataFrame(mc, index=output_top.index, columns=output_top.columns)

        # make a bedfile from the selected regions
        df = pd.DataFrame(output_mc, columns=['location'])
        df['location'] = df.index
        df['gene'] = df['location'].str.split(':').str[0]
        df['start'] = df['location'].str.split(':').str[1].str.split('-').str[0]
        df['end'] = df['location'].str.split(':').str[1].str.split('-').str[1]

        df = df.drop(['location'], axis=1)
        df = df.reset_index(drop=True)

        df.to_csv(r'/mbshome/mkolmus/danRer11_dif_top10000.bed', index = False, header=False, sep="\t")

        #export file as csv for gimme_maelstrom
        output_mc.to_csv(str(output), index = True, index_label = str, header=True, sep="\t") 

rule heatmap:
    input:
        bed_file = rules.transform_log_top_mc.output.output_bed,
        bam_files = get_inputfiles_bam
    params:
        clusters = config['clusters']
    output:
        config['out_dir'] + "heatmaps/{species}"
    shell:
        "fluff heatmap -f {input.bed_file} -d {input.bam_files} -o {output} -C k -k {params.clusters} -g -M Euclidean -c blue --no-colorbar"

rule gimme_maelstrom:
    input:
        table = rules.transform_log_top_mc.output,
        genome = config['in_dir'] + "genomes/{species}.fa"
    output:
        directory(config['out_dir'] + "{species}_maelstrom.out/")
    threads:
        48
    shell:
        "gimme maelstrom {input.table} {input.genome} {output} -N {threads}"
