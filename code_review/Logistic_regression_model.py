output = pd.DataFrame()

# these are all variables that can be tweaked to see which ones generate the best performing model:
param_grid = {
    "DEG_log2FC_cutoff": [0],
    "DEG_pvalue": [0.05],
    "H3K27ac_var_cutoff": [0],
    "H3K27ac_cutoff": [0],#[2.3],
    "H3K4me3_cutoff": [0],#[2.7],
    "H3K4me3_var_cutoff": [0.0],
    "distance_prom_regions": [20000],
    "min_counts":[100],
}

grid = ParameterGrid(param_grid)

# loop over all potential variable settings and run the code for all,
#append the output of each itteration to the 'output' list to compare how
#well different versions of the model performed.
for params in grid:

    # lets first find all differential genes based on the parameters of the model, using a
    # deseq2 result table and filtering on p-vallue
    DEG_pd = pd.read_table(
        DEG_file,
        sep=",",
        index_col=0,
        dtype={
            "baseMean": float,
            "log2FoldChange": float,
            "lfcSE": float,
            "stat": float,
            "pvalue": float,
            "padj": float,
            "gene": str,
            "total_counts": float,
        },
    )
    
    #filter away genes with to few counts
    DEG_pd = DEG_pd[DEG_pd.total_counts > params['min_counts']]

    # make a better annotation column, for now I give all DEGS a 1, however I sometimes play
    #arround with the type of prediction, therefore this 'overengeneerded' solution:
    conditions = [(
            (DEG_pd.log2FoldChange > params["DEG_log2FC_cutoff"])
            & (DEG_pd.padj < params["DEG_pvalue"]) 
        ),
        ((DEG_pd.log2FoldChange < -params["DEG_log2FC_cutoff"])
            & (DEG_pd.padj < params["DEG_pvalue"] )
        ),
    ]
    choices = [
        1,
        1,
    ]
    DEG_pd["gene_annotation"] = np.select(conditions, choices, default=0)

    # find the variable regions for the promoter and enhancer data, potentially I can 
    # set a threshold on how high a cis-regulatory element should have a signal
    # in order to be included in the analysis. 
    # however, this was adding pretty little info so right now both cutoffs are at 0.

    H3K27ac_regions = find_active_and_variable_regions(
        H3K27ac_norm,
        params["H3K27ac_cutoff"],
        params["H3K27ac_var_cutoff"],
        "KC_enh",
        "LSC_enh",
    )

    H3K27ac_var_regions = H3K27ac_regions[
        (
            (H3K27ac_regions.loc[:, "variable"] == True)
            & (H3K27ac_regions.loc[:, "active"] == True)
        )
    ]
    make_bedfile_from_column(
        H3K27ac_var_regions, "ensmbl_loc", f"{output_dir}/variable_H3K27ac.bed"
    )
    ##H3K4me3
    H3K4me3_regions = find_active_and_variable_regions(
        H3K4me3_norm,
        params["H3K4me3_cutoff"],
        params["H3K4me3_var_cutoff"],
        "KC_prom",
        "LSC_prom",
    )
    H3K4me3_var_regions = H3K4me3_regions[
        (
            (H3K4me3_regions.loc[:, "variable"] == True)
            & (H3K4me3_regions.loc[:, "active"] == True)
        )
    ]
    H3K4me3_var_regions[["chrom", "start", "end"]] = H3K4me3_var_regions[
        "ensmbl_loc"
    ].str.split("[:-]", expand=True)

    # Get the summit location of each H3K4me3 cis-regulatory region
    H3K4me3_var_regions["summit_start"] = (
        H3K4me3_var_regions["start"].astype(int) + H3K4me3_window / 2
    )
    H3K4me3_var_regions["summit_end"] = H3K4me3_var_regions["summit_start"] + 1
    H3K4me3_var_regions["summit_end"] = H3K4me3_var_regions["summit_end"].astype(int)
    H3K4me3_var_regions["summit_start"] = H3K4me3_var_regions["summit_start"].astype(
        int
    )

    H3K4me3_var_regions
    H3K4me3_var_regions["ensmbl_loc_summit"] = (
        H3K4me3_var_regions["chrom"]
        + ":"
        + H3K4me3_var_regions["summit_start"].astype(str)
        + "-"
        + H3K4me3_var_regions["summit_end"].astype(str)
    )
    make_bedfile_from_column(
        H3K4me3_var_regions, "ensmbl_loc_summit", f"{output_dir}/variable_H3K4me3.bed"
    )

    # link all TSS's to the variable bedfile locations. 
    # the h3K4me3 regions are maped to the closest TSS 
    # while the enhancers are mapped to all TSS regions within a 200kb window 
    
    TSS_2_H3K27ac = TSS_window_to_region(
        genome_path_gtf_small, f"{output_dir}/variable_H3K27ac.bed", f"100000"
    )
    TSS_2_H3K4me3 = TSS_to_region(
        genome_path_gtf_small, f"{output_dir}/variable_H3K4me3.bed", f"-k 1 -d"
    )

    # load the output bedfile with the included distance between each TSS and their closest variable H3K4me3 
    # histone modification,add the variable histone intensity score

    TSS_2_H3K27ac["loc"] = (
        TSS_2_H3K27ac["Chrom"]
        + ":"
        + TSS_2_H3K27ac["ChromStart"].astype(str)
        + "-"
        + TSS_2_H3K27ac["ChromEnd"].astype(str)
    )

    TSS_2_H3K27ac = TSS_2_H3K27ac.merge(
        H3K27ac_var_regions.iloc[:, [0, 1]],
        how="left",
        left_on="loc",
        right_index=True,
    )

    TSS_2_H3K4me3["loc"] = (
        TSS_2_H3K4me3["Chrom"]
        + ":"
        + TSS_2_H3K4me3["ChromStart"].astype(str)
        + "-"
        + TSS_2_H3K4me3["ChromEnd"].astype(str)
    )

    TSS_2_H3K4me3 = TSS_2_H3K4me3.merge(
        H3K4me3_var_regions.iloc[:, [0, 1, 12]],
        how="left",
        left_on="loc",
        right_on="ensmbl_loc_summit",
    )

    # calculate the intensity of each histone mark. For the enhancers use a distance weight matrix (from ANANSE).
    # For the H3K27me3, we will use a window acros the TSS, while for the
    # promoter marks we will map each TSS to the closest (variable) H3K4me3-ATAC region within 20kb.

    ##H3K27ac
    weighted_TSS_2_H3K27ac = distance_weight_region_average(TSS_2_H3K27ac, weight_dict)
    weighted_TSS_2_H3K27ac["gene_name"] = weighted_TSS_2_H3K27ac.index
    weighted_TSS_2_H3K27ac = weighted_TSS_2_H3K27ac.rename(
        columns={
            0: "KC",
            1: "LSC",
            2: "n_enh",
            "mean_int": "mean_int_ac",
            "abs_FC": "abs_FC_ac",
            "FC": "FC_ac",
        }
    )
    weighted_TSS_2_H3K27ac = weighted_TSS_2_H3K27ac[
        ["gene_name", "mean_int_ac", "FC_ac", "abs_FC_ac", "n_enh"]
    ]

    ## calculate the mean H3K4me3 signal and the (absolute) FC between KC and LSC
    TSS_2_H3K4me3["mean_int_prom"] = np.mean(
        TSS_2_H3K4me3.loc[:, ["KC", "LSC"]], axis=1
    )
    TSS_2_H3K4me3["FC_prom"] = np.log2(
        TSS_2_H3K4me3["KC"].astype(float) / TSS_2_H3K4me3["LSC"]
    ).astype(float)
    TSS_2_H3K4me3["FC_abs_prom"] = abs(TSS_2_H3K4me3["FC_prom"])

    averaged_TSS_2_H3K4me3 = TSS_2_H3K4me3[
        ["gene_name", "mean_int_prom", "FC_prom", "FC_abs_prom", "gene_region_dist"]
    ]
    averaged_TSS_2_H3K4me3 = averaged_TSS_2_H3K4me3.drop_duplicates(subset="gene_name")
    averaged_TSS_2_H3K4me3 = averaged_TSS_2_H3K4me3[
        averaged_TSS_2_H3K4me3["gene_region_dist"] < params["distance_prom_regions"]
    ]

    TSS_2_all = weighted_TSS_2_H3K27ac.merge(
        averaged_TSS_2_H3K4me3, how="outer", on="gene_name"
    ).merge(averaged_TSS_2_H3K27me3, how="outer", left_on="gene_name", right_index=True)

    TSS_2_all = TSS_2_all.merge(
        DEG_pd[["log2FoldChange", "padj", "gene_annotation", "baseMean"]],
        how="left",
        left_on="gene_name",
        right_index=True,
    )

    TSS_2_all = TSS_2_all[TSS_2_all["gene_annotation"].notna()]
    TSS_2_all = TSS_2_all.replace([np.inf, -np.inf], np.nan)


    # lets make the X and y df for the logistic regression model:
    X = TSS_2_all.loc[
        :,
        [
            "gene_name",
            "mean_int_ac",
            "abs_FC_ac",
            "n_enh",
            "mean_int_prom",
            "FC_abs_prom",
            "gene_region_dist",
            "mean_int_repr",
            "FC_abs_repr",
        ],
    ]
    X = X.astype(
        {
            "gene_name": str,
            "mean_int_ac": np.float32,
            "abs_FC_ac": np.float32,
            "n_enh": np.float32,
            "mean_int_prom": np.float32,
            "FC_abs_prom": np.float32,
            "gene_region_dist": np.float32,
            "mean_int_repr": np.float32,
            "FC_abs_repr": np.float32,
        }
    )
    y = TSS_2_all.loc[:, ["gene_annotation", "gene_name"]]

    # next we will split the dataset in test and train data, using all genes on chrom 1 as a test dataset
    chr1_gene_names = gtf_df[gtf_df["Chrom_TSS"] == "1"]["gene"]

    X_test = X[X["gene_name"].isin(chr1_gene_names)].drop("gene_name", axis=1)
    y_test = y[y["gene_name"].isin(chr1_gene_names)].drop("gene_name", axis=1)
    X_train = X[-X["gene_name"].isin(chr1_gene_names)].drop("gene_name", axis=1)
    y_train = y[-y["gene_name"].isin(chr1_gene_names)].drop("gene_name", axis=1)

    X = X.drop("gene_name", axis=1)
    y = y.drop("gene_name", axis=1)

    # since it is an unbalanced dataste (way more genes are non-diff than diff), we will calculate the ratio and
    # feed this to the model so it gives a higher weight to guessing the diff genes right.

    weights = {
        1: TSS_2_all["gene_annotation"].value_counts(normalize=True)[0],
        0: TSS_2_all["gene_annotation"].value_counts(normalize=True)[1],
    }

    # Create a logistic regression pipeline,
    # first we impute values for the NaN vallues present in the DF, then wer
    #perform scaling (mostly for the n_enh, and distance to promoter values)
    #while finally running logistic regression
    steps = [
        (
            "imputation",
            SimpleImputer(missing_values=np.nan, strategy='most_frequent', fill_value=0),
        ),
        ("scaler", StandardScaler()),
        (
            "model",
            LogisticRegression(
                solver="liblinear", class_weight=weights, penalty="l1", C=1,
            ),
        ),
    ]
    pipeline = Pipeline(steps)
    
    #train the entire pipeline based on all genes except the ones on chr1
    pipeline.fit(
        X_train, y_train,
    )
    
    #calculate the ROC curve with 3 different subset from all the data
    cv_auc = cross_val_score(pipeline, X, y, cv=3, scoring="roc_auc")

    #predict if the genes on chr1 are differential or not
    y_pred = pipeline.predict(X_test)
    y_pred_prob = pipeline.predict_proba(X_test)[:, 1]
    
    #how precise was the model if we compare 
    av_PC = average_precision_score(y_test, y_pred_prob)

    #save all the QC metrics of the settings with this model.
    result = {
        "cutoff_vallues": str(params),
        "accuracy": str(metrics.accuracy_score(y_test, y_pred)),
        "cv_AUC": np.mean(cv_auc),
        "average_precision_score": av_PC,
    }
    print(result)
    output = output.append(result, ignore_index=True)
