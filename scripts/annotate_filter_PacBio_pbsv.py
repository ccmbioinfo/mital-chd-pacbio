import argparse
import os
import pandas as pd
from pybedtools import BedTool
from sigfig import round


def create_chd_gene_dict(chd):
    chd_entrez_dict = {}
    for index, row in chd.iterrows():
        chd_entrez_dict[row["EntrezID"]] = row["Class"]
    return chd_entrez_dict


def chd_gene(sv_entrez_id, chd_entrez_dict):
    try:
        sv_entrez_id = sv_entrez_id.replace(",", "|")
        sv_entrez_id = sv_entrez_id.split("|")
        class_list = []
        for gene in sv_entrez_id:
            gene = int(gene)
            if gene in chd_entrez_dict:
                c = chd_entrez_dict[gene]
                class_list.append(c)
            else:
                class_list.append("NA")
        return ("|").join(class_list)
    except AttributeError:
        # SV does not overlap any genes
        return "NA"


def sv_to_bed(sv):
    sv = sv[["CHROM", "POS", "END", "SVTYPE", "ID"]]
    bed = BedTool.from_dataframe(sv)
    return bed


def annotate_pop_svs(sv, pop_svs, cols):
    # intersect annotsv and population SV bed file
    sv_bed = sv_to_bed(sv)
    pop_svs = pd.read_csv(pop_svs, sep="\t")
    pop_svs_bed = BedTool.from_dataframe(pop_svs)
    intersect = sv_bed.intersect(
        pop_svs_bed, wa=True, wb=True, F=0.5, f=0.5
    ).to_dataframe()
    intersect.columns = [
        "CHROM",
        "POS",
        "END",
        "SVTYPE",
        "ID",
        "CHROM_pop",
        "POS_pop",
        "END_pop",
        "SVTYPE_pop",
    ] + cols
    # popSV and sample SV must be same type
    intersect = intersect[intersect["SVTYPE"] == intersect["SVTYPE_pop"]]
    pop_name = cols[0].split("_")[0]
    # e.g. make a column with SV details, e.g DEL:1:25266309-25324509
    intersect[f"{pop_name}_SV"] = intersect[
        ["CHROM_pop", "POS_pop", "END_pop", "SVTYPE_pop"]
    ].apply(lambda x: f"{x[3]}:{x[0]}:{x[1]}-{x[2]}", axis=1)
    cols.append(f"{pop_name}_SV")
    intersect = intersect[["CHROM", "POS", "END", "SVTYPE", "ID"] + cols]
    # round AFs
    try:
        AF_col = [col for col in cols if "AF" in col][0]
        intersect[AF_col] = [
            round(float(af), sigfigs=3) for af in intersect[AF_col].values
        ]
    except IndexError:
        pass
    # group by SV, joining annotation columns
    intersect = intersect.astype(str)
    intersect = (
        intersect.groupby(["CHROM", "POS", "END", "SVTYPE", "ID"])[cols]
        .agg({col: "; ".join for col in cols})
        .reset_index()
    )
    # get max allele frequency
    try:
        AF_col = [col for col in cols if "AF" in col][0]
        intersect[f"{pop_name}_maxAF"] = intersect[AF_col].apply(
            lambda x: max([float(af) for af in x.split("; ")])
        )
        cols.append(f"{pop_name}_maxAF")
    # get max allele counts for C4R
    except IndexError:
        count_cols = ["C4R_AC", "seen_in_C4R_count"]
        for col in count_cols:
            intersect[f"{col}_max"] = intersect[col].apply(
                lambda x: max([int(ac) for ac in x.split("; ")])
            )
        cols.append(f"{col}_max")

    # merge population AF dataframe with SV df
    sv_pop_svs = pd.merge(
        sv,
        intersect,
        how="left",
        on=["CHROM", "POS", "END", "SVTYPE", "ID"],
    ).fillna(value={col: 0 for col in cols})
    return sv_pop_svs


def main(sv, cmh, chd, output):
    sv_filename =os.path.join(output,os.path.basename(sv).replace(".tsv", ""))
    sv = pd.read_csv(sv, sep="\t")
    chd = pd.read_csv(args.chd, sep="\t")
    chd_entrez_dict = create_chd_gene_dict(chd)

    # add CHD genes
    sv["CHD_class"] = sv.loc[:, "gene_egID"].apply(
        lambda row: chd_gene(row, chd_entrez_dict)
    )

    # filter out non-CHD genes
    sv_chd = sv[sv["CHD_class"].str.contains("Candidate|Modifier|Syndromic", regex = True)]

    # filter out SVs < 50 bp
    sv_over_50 = sv_chd[sv_chd["SVLEN"] >= 50]

    # remove 'chr' prefix for compatibility with CMH bed file
    sv_over_50["CHROM"] = sv_over_50.loc[:, "CHROM"].apply(lambda x: x.replace("chr", ""))
    sv_over_50 = sv_over_50.astype(str)
    cmh_cols = ["cmh_AF", "cmh_AC"]

    # annotate svs with CMH allele frequencies and allele counts
    sv_cmh_over_50 = annotate_pop_svs(sv_over_50, cmh, cmh_cols)

    # filter for rare svs
    sv_cmh_over_50 = sv_cmh_over_50.astype({"pacBioPercFreq_90percRecipOverlap": float, "cmh_maxAF": float})
    # PacBio is allele frequency, CMH is allele fraction
    # don't filter on PacBio for now, don't know how many individuals are in TCAG's cohort
    # sv_cmh_over_50_rare = sv_cmh_over_50[
    #     (sv_cmh_over_50["pacBioPercFreq_90percRecipOverlap"] < 1) & (sv_cmh_over_50["cmh_maxAF"] < 0.01)
    # ]
    sv_cmh_over_50_rare = sv_cmh_over_50[ sv_cmh_over_50["cmh_maxAF"] < 0.01]

    # output csvs
    sv_cmh_over_50_rare.to_csv(f"{sv_filename}_CHD_over_50_rare.csv", index=False)
    print("DONE")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates a rare structural variant report using TCAG-annotated pbsv calls generated from PacBio HiFi WGS"
    )
    parser.add_argument("-sv", type=str, help="TCAG pbsv tsv file", required=True)
    parser.add_argument(
        "-cmh",
        type=str,
        help="Children's Mercy Hospital pbsv allele counts and frequencies",
        required=True,
    )
    parser.add_argument(
        "-chd",
        help="Tab delimited file containing gene names and HPO terms",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-output",
        help="Output directory to write annotated reports",
        type=str,
        required=True,
    )

args = parser.parse_args()
main(args.sv, args.cmh, args.chd, args.output)

