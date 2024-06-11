# %%
import pandas as pd

# %%
fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
hsdf = pd.read_csv(fname, index_col=0)
fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_pangenome.csv"
hspg = pd.read_csv(fname, index_col=0)

# %%
x = pd.DataFrame(
    [
        [18, "CIRMBUYJFK_f__CWCCKOQCWZ_r"],
        [20, "HEDBBNUKLU_r__VSCDQDGRHO_r"],
        [39, "ETMEUQAZWU_f__QSPABZAPOJ_f"],
        [36, "SFNIQHXIST_r__TJOBMLQRFA_r"],
        [12, "BWEZXGGFBK_r__MVMOFPVELT_r"],
        [21, "IIFSSHHHUK_f__RJJWLWHZAS_f"],
        [19, "TLVFRBMGBC_r__YUMHUOWTXQ_f"],
        [31, "ATFHQYFPNW_f__GUDQDOMJFJ_f"],
        [24, "GPKQYOCEJI_r__NKVSUZGURN_f"],
        [10, "EQDOECTTHL_f__FSKJDTCZAX_r"],
        [40, "ATPWUNKKID_f__KKPYPKGMXA_f"],
        [35, "LGOKMQPFVQ_f__XNYZXWCUST_r"],
        [17, "KAIPBNCIHR_r__WXCHSHHCDT_f"],
        [11, "RKAOKULCFF_f__VFXLTFSKTV_r"],
        [37, "RYYAQMEJGY_f__XKJZBXCDPZ_f"],
        [37, "RYYAQMEJGY_r__ZTHKZYHPIX_f"],
        [27, "IHKFSQQUKE_r__KPBYGJHRZJ_f"],
        [28, "IHKFSQQUKE_r__KPBYGJHRZJ_f"],
        [29, "IHKFSQQUKE_r__KPBYGJHRZJ_f"],
        [38, "JVNRLCFAVD_f__PLTCZQCVRD_r"],
        [5, "XXVMWZCEKI_r__YUOECYBHUS_r"],
    ],
    columns=["hsn", "jn"],
)
x = x.join(hsdf["n_categories"], on="jn")
x = x.join(hspg["pangenome_len"], on="jn")
x["pgl"] = (x["pangenome_len"] / 1000).round(0)
x
# %%
for idx, row in x.iterrows():
    print(
        f"[{row.hsn}], [{row.jn}], [{row.pgl} kbp], [{row.n_categories}], [], [], [],"
    )

# %%
