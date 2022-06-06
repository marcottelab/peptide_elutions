
############ Artificial trypsin digest

#June 2017 Uniprot proteome download
#python scripts/trypsin.py --input proteomes/arath/uniprot-proteome%3AUP000006548.fasta  --output  proteomes/arath/uniprot-proteome%3AUP000006548_peptides.csv  --miss 2
#python scripts/trypsin.py --input proteomes/braol/uniprot-proteome%3AUP000032141.fasta  --output  proteomes/braol/uniprot-proteome%3AUP000032141_peptides.csv  --miss 2
#python scripts/trypsin.py --input proteomes/chlre/uniprot-proteome%3AUP000006906.fasta  --output  proteomes/chlre/uniprot-proteome%3AUP000006906_peptides.csv --miss 2
#python scripts/trypsin.py --input proteomes/orysj/uniprot-proteome%3AUP000059680.fasta  --output  proteomes/orysj/uniprot-proteome%3AUP000059680_peptides.csv --miss 2
#python scripts/trypsin.py --input proteomes/phatc/uniprot-proteome%3AUP000000759.fasta  --output  proteomes/phatc/uniprot-proteome%3AUP000000759_peptides.csv --miss 2
#python scripts/trypsin.py --input proteomes/selml/uniprot-proteome%3AUP000001514.fasta  --output  proteomes/selml/uniprot-proteome%3AUP000001514_peptides.csv --miss 2
#python scripts/trypsin.py --input proteomes/traes/uniprot-proteome%3AUP000019116.fasta  --output  proteomes/traes/uniprot-proteome%3AUP000019116_peptides.csv --miss 2
#python scripts/trypsin.py --input proteomes/lytva/Lvar2_2.fasta --output proteomes/lytva/Lvar2_2_peptides.csv --miss 2 
#python scripts/trypsin.py --input proteomes/chqui/Cquinoa_392_v1.0.protein.fasta --output proteomes/chqui/Cquinoa_392_v1.0_peptides.csv --miss 2 

python scripts/trypsin.py --input proteomes/mouse/uniprot-proteome%3AUP000000589.fasta  --output  proteomes/mouse/uniprot-proteome%3AUP000000589_peptides.csv --miss 2
python scripts/trypsin.py --input proteomes/human/uniprot-proteome%3AUP000005640.fasta  --output  proteomes/human/uniprot-proteome%3AUP000005640_peptides.csv --miss 2
python scripts/trypsin.py --input proteomes/caeel/uniprot-proteome%3AUP000001940.fasta  --output  proteomes/caeel/uniprot-proteome%3AUP000001940_peptides.csv --miss 2
python scripts/trypsin.py --input proteomes/dicdi/uniprot-proteome%3AUP000002195.fasta  --output  proteomes/dicdi/uniprot-proteome%3AUP000002195_peptides.csv --miss 2
python scripts/trypsin.py --input proteomes/drome/uniprot-proteome%3AUP000000803.fasta  --output  proteomes/drome/uniprot-proteome%3AUP000000803_peptides.csv --miss 2

python scripts/trypsin.py --input proteomes/yeast/uniprot-proteome%3AUP000002311.fasta  --output  proteomes/yeast/uniprot-proteome%3AUP000002311_peptides.csv --miss 2


python scripts/trypsin.py --input proteomes/xenla/XENLA_JGIv18p_prot_annot_HS.fasta  --output  proteomes/xenla/XENLA_JGIv18p_prot_annot_HS_peptides.csv --miss 2



python scripts/trypsin.py --input proteomes/sollc/uniprot-proteome%3AUP000004994.fasta  --output  proteomes/sollc/uniprot-proteome%3AUP000004994_peptides.csv --miss 2
python scripts/trypsin.py --input proteomes/tetts/uniprot-proteome%3AUP000009168.fasta  --output  proteomes/tetts/uniprot-proteome%3AUP000009168_peptides.csv --miss 2

python scripts/trypsin.py --input_file proteomes/strpu/working_proteome/uniprot-proteome%3AUP000007110.fasta  --output_file  proteomes/strpu/working_proteome/uniprot-proteome%3AUP000007110_peptides.csv --miss 2

python scripts/trypsin.py --input_file proteomes/nemve/working_proteome/uniprot-proteome%3AUP000001593.fasta  --output_file  proteomes/nemve/working_proteome/uniprot-proteome%3AUP000001593_peptides.csv --miss 2
python scripts/trypsin.py --input_file proteomes/xenla/working_proteome/uniprot-proteome%3AUP000186698.fasta  --output_file  proteomes/xenla/working_proteome/uniprot-proteome%3AUP000186698_peptides.csv --miss 2


########### DIAMOND lookup
#bash record_diamond.sh




########### Define grouping for peptide-protein lookup
python scripts/define_grouping.py --spec arath --grouping_type protein  --peptides proteomes/arath/uniprot-proteome%3AUP000006548_peptides.csv --output_dir proteomes/arath/
python scripts/define_grouping.py --spec arath --grouping_type euNOG --grouping eggnog_mapper/arath_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/arath/uniprot-proteome%3AUP000006548_peptides.csv --output_dir proteomes/arath/


python scripts/define_grouping.py --spec braol --grouping_type protein  --peptides proteomes/braol/uniprot-proteome%3AUP000032141_peptides.csv --output_dir proteomes/braol/
python scripts/define_grouping.py --spec braol --grouping_type euNOG --grouping eggnog_mapper/braol_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/braol/uniprot-proteome%3AUP000032141_peptides.csv --output_dir proteomes/braol/


python scripts/define_grouping.py --spec chlre --grouping_type protein  --peptides proteomes/chlre/uniprot-proteome%3AUP000006906_peptides.csv --output_dir proteomes/chlre/
python scripts/define_grouping.py --spec chlre --grouping_type euNOG --grouping eggnog_mapper/chlre_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/chlre/uniprot-proteome%3AUP000006906_peptides.csv --output_dir proteomes/chlre/


python scripts/define_grouping.py --spec orysj --grouping_type protein  --peptides proteomes/orysj/uniprot-proteome%3AUP000059680_peptides.csv --output_dir proteomes/orysj/
python scripts/define_grouping.py --spec orysj --grouping_type euNOG --grouping eggnog_mapper/orysj_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/orysj/uniprot-proteome%3AUP000059680_peptides.csv --output_dir proteomes/orysj/


python scripts/define_grouping.py --spec phatc --grouping_type protein  --peptides proteomes/phatc/uniprot-proteome%3AUP000000759_peptides.csv --output_dir proteomes/phatc/
python scripts/define_grouping.py --spec phatc --grouping_type euNOG --grouping eggnog_mapper/phatc_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/phatc/uniprot-proteome%3AUP000000759_peptides.csv --output_dir proteomes/phatc/


python scripts/define_grouping.py --spec selml --grouping_type protein  --peptides proteomes/selml/uniprot-proteome%3AUP000001514_peptides.csv --output_dir proteomes/selml/
python scripts/define_grouping.py --spec selml --grouping_type euNOG --grouping eggnog_mapper/selml_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/selml/uniprot-proteome%3AUP000001514_peptides.csv --output_dir proteomes/selml/


python scripts/define_grouping.py --spec traes --grouping_type protein  --peptides proteomes/traes/uniprot-proteome%3AUP000019116_peptides.csv --output_dir proteomes/traes/
python scripts/define_grouping.py --spec traes --grouping_type euNOG --grouping eggnog_mapper/traes_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/traes/uniprot-proteome%3AUP000019116_peptides.csv --output_dir proteomes/traes/

python scripts/define_grouping.py --spec sollc --grouping_type protein  --peptides proteomes/sollc/uniprot-proteome%3AUP000004994_peptides.csv --output_dir proteomes/sollc/
python scripts/define_grouping.py --spec sollc --grouping_type euNOG --grouping eggnog_mapper/sollc_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/sollc/uniprot-proteome%3AUP000004994_peptides.csv --output_dir proteomes/sollc/


python scripts/define_grouping.py --spec tetts --grouping_type protein  --peptides proteomes/tetts/uniprot-proteome%3AUP000009168_peptides.csv --output_dir proteomes/tetts/
python scripts/define_grouping.py --spec tetts --grouping_type euNOG --grouping eggnog_mapper/tetts_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/tetts/uniprot-proteome%3AUP000009168_peptides.csv --output_dir proteomes/tetts/


python scripts/define_grouping.py --spec mouse --grouping_type protein  --peptides proteomes/mouse/uniprot-proteome%3AUP000000589_peptides.csv --output_dir proteomes/mouse/

python scripts/define_grouping.py --spec mouse --grouping_type euNOG --grouping eggnog_mapper/mouse_hmmer_euk.mapping  --peptides proteomes/mouse/uniprot-proteome%3AUP000000589_peptides.csv --output_dir proteomes/mouse/


python scripts/define_grouping.py --spec human --grouping_type protein  --peptides proteomes/human/uniprot-proteome_human_reviewed_peptides.csv --output_dir proteomes/human/
python scripts/define_grouping.py --spec human --grouping_type euNOG --grouping eggnog_mapper/human_hmmer_euk.mapping  --peptides proteomes/human/uniprot-proteome_human_reviewed_peptides.csv --output_dir proteomes/human/


python scripts/define_grouping.py --spec caeel --grouping_type protein  --peptides proteomes/caeel/uniprot-proteome%3AUP000001940_peptides.csv --output_dir proteomes/caeel/
python scripts/define_grouping.py --spec caeel --grouping_type euNOG --grouping eggnog_mapper/caeel_hmmer_euk.mapping  --peptides proteomes/caeel/uniprot-proteome%3AUP000001940_peptides.csv --output_dir proteomes/caeel/

python scripts/define_grouping.py --spec dicdi --grouping_type protein  --peptides proteomes/dicdi/uniprot-proteome%3AUP000002195_peptides.csv --output_dir proteomes/dicdi/
python scripts/define_grouping.py --spec dicdi --grouping_type euNOG --grouping eggnog_mapper/dicdi_hmmer_euk.mapping  --peptides proteomes/dicdi/uniprot-proteome%3AUP000002195_peptides.csv --output_dir proteomes/dicdi/

python scripts/define_grouping.py --spec drome --grouping_type protein  --peptides proteomes/drome/uniprot-proteome%3AUP000000803_peptides.csv --output_dir proteomes/drome/
python scripts/define_grouping.py --spec drome --grouping_type euNOG --grouping eggnog_mapper/drome_hmmer_euk.mapping  --peptides proteomes/drome/uniprot-proteome%3AUP000000803_peptides.csv --output_dir proteomes/drome/


python scripts/define_grouping.py --spec yeast --grouping_type protein  --peptides proteomes/yeast/uniprot-proteome%3AUP000002311_peptides.csv --output_dir proteomes/yeast/
python scripts/define_grouping.py --spec yeast --grouping_type euNOG --grouping eggnog_mapper/yeast_hmmer_euk.mapping  --peptides proteomes/yeast/uniprot-proteome%3AUP000002311_peptides.csv --output_dir proteomes/yeast/

python scripts/define_grouping.py --spec xenla --grouping_type protein  --peptides proteomes/xenla/XENLA_JGIv18p_prot_annot_HS_peptides.csv --output_dir proteomes/xenla/
python scripts/define_grouping.py --spec xenla --grouping_type euNOG --grouping eggnog_mapper/xenla_hmmer_euk.mapping  --peptides proteomes/xenla/XENLA_JGIv18p_prot_annot_HS_peptides.csv --output_dir proteomes/xenla/


python scripts/define_grouping.py --spec strpu --grouping_type protein  --peptides proteomes/strpu/working_proteome/uniprot-proteome%3AUP000007110_peptides.csv --output_dir proteomes/strpu/
python scripts/define_grouping.py --spec strpu --grouping_type euNOG --grouping eggnog_mapper/strpu_hmmer_euk.mapping  --peptides proteomes/strpu/uniprot-proteome%3AUP000007110_peptides.csv --output_dir proteomes/strpu/

python scripts/define_grouping.py --spec nemve --grouping_type protein  --peptides proteomes/nemve/working_proteome/uniprot-proteome%3AUP000001593_peptides.csv --output_dir proteomes/nemve/
python scripts/define_grouping.py --spec nemve --grouping_type euNOG --grouping eggnog_mapper/nemve_hmmer_euk.mapping  --peptides proteomes/nemve/uniprot-proteome%3AUP000001593_peptides.csv --output_dir proteomes/nemve/

python scripts/define_grouping.py --spec xenla --grouping_type protein  --peptides proteomes/xenla/working_proteome/uniprot-proteome%3AUP000186698_peptides.csv --output_dir proteomes/xenla/
python scripts/define_grouping.py --spec xenla --grouping_type euNOG --grouping eggnog_mapper/xenla_hmmer_euk.mapping  --peptides proteomes/xenla/uniprot-proteome%3AUP000186698_peptides.csv --output_dir proteomes/xenla/


python scripts/define_grouping.py --spec arath --grouping_type euNOG --grouping eggnog_mapper/arath_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/arath/uniprot-proteome%3AUP000006548_peptides.csv --output_dir proteomes/arath/

python scripts/define_grouping.py --spec braol --grouping_type euNOG --grouping eggnog_mapper/braol_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/braol/uniprot-proteome%3AUP000032141_peptides.csv --output_dir proteomes/braol/

python scripts/define_grouping.py --spec chlre --grouping_type euNOG --grouping eggnog_mapper/chlre_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/chlre/uniprot-proteome%3AUP000006906_peptides.csv --output_dir proteomes/chlre/

python scripts/define_grouping.py --spec orysj --grouping_type euNOG --grouping eggnog_mapper/orysj_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/orysj/uniprot-proteome%3AUP000059680_peptides.csv --output_dir proteomes/orysj/

python scripts/define_grouping.py --spec selml --grouping_type euNOG --grouping eggnog_mapper/selml_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/selml/uniprot-proteome%3AUP000001514_peptides.csv --output_dir proteomes/selml/

python scripts/define_grouping.py --spec traes --grouping_type euNOG --grouping eggnog_mapper/traes_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/traes/uniprot-proteome%3AUP000019116_peptides.csv --output_dir proteomes/traes/

python scripts/define_grouping.py --spec sollc --grouping_type euNOG --grouping eggnog_mapper/sollc_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/sollc/uniprot-proteome%3AUP000004994_peptides.csv --output_dir proteomes/sollc/

python scripts/define_grouping.py --spec chqui --grouping_type euNOG  --peptides proteomes/chqui/Cquinoa_392_v1.0_peptides.csv --output_dir proteomes/chqui/


:wq

############Convert from pep_count to matrix
######All at once
#for f in pep_counts_folder/??_*; do bash scripts/consolidate_MSblender_output.sh $f ${f##*/} elutions; done
######One at a time
#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_Ptricornutum_IEX_20170328/output OP_Ptricornutum_IEX_20170328 elutions/
#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_Ptricornutum_SEC_20170403/output OP_Ptricornutum_SEC_20170403 elutions/

#bash scripts/consolidate_MSblender_output.sh /MS/processed/Lumos/Ophelia/Euglena_SEC_1-2017/output Euglena_SEC_1-2017 elutions/
#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_Arabidopsis_sproutsSEC_02102017/output OP_Arabidopsis_sproutsSEC_02102017 elutions/
#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/RiceLeaf_SEC_04292016/output RiceLeaf_SEC_04292016 elutions/
#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_WheatGrassSEC_20170202/output OP_WheatGrassSEC_20170202 elutions/
#bash scripts/consolidate_MSblender_output.sh /MS/processed/Lumos/Arabidopsis_sprouts_WWC_2017/msblender Arabidopsis_sprouts_WWC_2017 elutions/ 
#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_SeaUrchinSperm_WWC_20160119/msblender OP_SeaUrchinSperm_WWC_20160119 elutions/ 
#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_Chlamy_Nuclei_SEC_20170428/output OP_Chlamy_Nuclei_SEC_20170428 elutions/ 
#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_Quinoa_IEX_20170513/output OP_Quinoa_IEX_20170513 elutions/ 

#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_Tetrahymena_bodies_IEX_20170619/output OP_Tetrahymena_bodies_IEX_20170619 elutions/ 

bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/Anna_HEK293T_IEX_control_011217/output Anna_HEK293T_IEX_control_011217 elutions/ 

bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/Anna_HEK293T_IEX_mixed_bed_control_042217/output Anna_HEK293T_IEX_mixed_bed_control_042217 elutions/ 

bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/Anna_HEK293T_IEX_mixed_bed_RNASE_050117/output Anna_HEK293T_IEX_mixed_bed_RNASE_050117 elutions/ 

bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/Anna_HEK293T_IEX_RNaseA_012517/output Anna_HEK293T_IEX_RNaseA_012517 elutions/ 

bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/Anna_HEK293T_SEC_Benzonase_121316/output Anna_HEK293T_SEC_Benzonase_121316 elutions/ 

bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/Anna_HEK293T_SEC_control_113016/output Anna_HEK293T_SEC_control_113016 elutions/ 

bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/Anna_HEK293T_SEC_RNaseA_120916/output Anna_HEK293T_SEC_RNaseA_120916 elutions/ 

bash script/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_Tetrahymena_cilia_SEC_20170708/output OP_Tetrahymena_cilia_SEC_20170708 elutions/


bash scripts/consolidate_MSblender_output.sh /MS/processed/Lumos/RBC_1086_WWC_ghosts_prep1_11232016/output RBC_1086_WWC_ghosts_prep1_11232016 elutions/ 

for exp in Xl_caps_1308 Xl_musc_op_22_1303 Xl_musc_op_26_1212 Xl_muscrna_op_24_1302
do 
    bash scripts/consolidate_MSblender_output.sh /MS/processed/metazoans/$exp/output $exp elutions/ 
done

for exp in Sp_60_121010 Sp_beadsALF_60_121014 Sp_beadsBNF_73_1201119 Sp_beadsL_60_121125 Sp_iex_60_121019 
do 
    bash scripts/consolidate_MSblender_output.sh /MS/processed/metazoans/$exp/output $exp elutions/ 
done

for exp in  Nv_hcw_110_1207 Nv_iex_70_121224 Nv_iex_70_121229 Nv_iex_82_130111 Nv_iex_91_130103 Nv_iex_96_121213  
do 
    bash scripts/consolidate_MSblender_output.sh /MS/processed/metazoans/$exp/output $exp elutions/ 
done




for exp in Ce_1104 Ce_6mg_1202 Ce_6mg_1203 Ce_beadsA_1206 Ce_beadsB_wan_60_1208 Ce_beadsC_wan_60_1209 Ce_beadsflow_1206 Ce_beadsL_1206 Ce_beadsN_wan_60_1208 Ce_beadsR_wan_60_1209 Ce_BNF_wan_60_1209 Ce_CRF_wan_60_1209 Ce_ups_1205 Ce_ups_1206 
do 
    bash scripts/consolidate_MSblender_output.sh /project/cmcwhite/split_profiles/pep_counts_folder/$exp $exp elutions/ 
done


for exp in Ec_quiestran_JDO_201307
do 
    bash scripts/consolidate_MSblender_output.sh /project/cmcwhite/split_profiles/pep_counts_folder/$exp $exp elutions/ 
done



for exp in Hs_helaC_ph_tcs269_2 Hs_helaC_ph_tcs269_1 Hs_helaC_mar_SGF Hs_hekN_1201 Hs_helaN_ph_tcs375_P2 Hs_hekN_ph_hcw120_1 Hs_helaN_ph_tcs375_P3 Hs_helaN_mar_ief3_10 Hs_helaN_mar_ief5_8 Hs_CB660_1105 Hs_helaN_mar_SGF Hs_G166_1104 Hs_helaN_ph_hcw120_1 Hs_G166_1105 Hs_helaN_ph_hcw120_2 Hs_hek_kinase95 Hs_helaN_ph_saf48 Hs_hekN_1108 Hs_helaN_ph_tcs375_P1 Hs_hekN_ph_hcw120_2 Hs_helaN_ph_tcs375_P4 Hs_hek_wan_membraneopioid72 Hs_helaN_wan_hcw120_3 Hs_helaC_mar_ief3_10 Hs_wan_aguhplc108_1204 Hs_helaC_mar_ief3_10_repeat Hs_wan_shuhplc94_1206 Hs_helaC_mar_ief5_8 Hs_helaC_ph_wax43 Hs_helaN_1003 Hs_helaN_1010 
do 
    bash scripts/consolidate_MSblender_output.sh /project/cmcwhite/split_profiles/pep_counts_folder/$exp $exp elutions/ 
done


for exp in yeast?


for exp in Mm_1106 Mm_1201
do 
    bash scripts/consolidate_MSblender_output.sh /project/cmcwhite/split_profiles/pep_counts_folder/$exp $exp elutions/ 
done

for exp in Dd_1110 Dd_1112 Dm_1107_Dem Dm_1107_DSL2 Dm_1201_Dem 
do 
    bash scripts/consolidate_MSblender_output.sh /project/cmcwhite/split_profiles/pep_counts_folder/$exp $exp elutions/ 
done



for exp in Xl_muscrna_op_24_1302 Xl_caps_1308 Xl_musc_op_26_1212 Xl_musc_op_22_1303 
do 
    bash scripts/consolidate_MSblender_output.sh /project/cmcwhite/split_profiles/pep_counts_folder/$exp $exp elutions/ 
done

#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_Tomato_leaf_IEX2_201270607/output OP_Tomato_leaf_IEX2_201270607 elutions/ 

#bash scripts/consolidate_MSblender_output.sh /MS/processed/Fusion_data/OP_Tetrahymena_cilia_IEX_20170609/output/ OP_Tetrahymena_cilia_IEX_20170609 elutions/ 



############Define peptide grouping for a species with a particular definition
###This takes a file with GroupID ProteinID columns and a peptide lookup with ProteinID and Peptide columns

#python scripts/define_grouping.py --spec euggr --grouping_type euNOG --grouping eggnog_output/euggr.Euglena_nr.euNOG_orthology.tab --peptides proteomes/euggr/Euglena_nr_peptides.csv --output_dir proteomes/euggr/

#python scripts/define_grouping.py --spec arath --grouping_type euNOG --grouping eggnog_output/arath.euNOG_orthology.tab --peptides proteomes/arath/uniprot-proteome%3AUP000006548_peptides.csv --output_dir proteomes/arath/
#python scripts/define_grouping.py --spec orysj --grouping_type euNOG --grouping eggnog_output/orysj.euNOG_orthology.tab --peptides proteomes/orysj/uniprot-proteome%3AUP000000763_peptides.csv --output_dir proteomes/orysj/
#python scripts/define_grouping.py --spec traes --grouping_type euNOG --grouping eggnog_output/traes.euNOG_orthology.tab --peptides proteomes/traes/uniprot-triticum%3AUP000019116_peptides.csv --output_dir proteomes/traes/
#python scripts/define_grouping.py --spec braol --grouping_type euNOG --grouping eggnog_output/braol.euNOG_orthology.tab --peptides proteomes/braol/uniprot-brassica+oleracea_peptides.csv --output_dir proteomes/braol/
#python scripts/define_grouping.py --spec selml --grouping_type euNOG --grouping eggnog_output/selml.euNOG_orthology.tab --peptides proteomes/selml/uniprot-proteome%3AUP000001514_peptides.csv --output_dir proteomes/selml/

#python scripts/define_grouping.py --spec chlre --grouping_type euNOG --grouping eggnog_output/chlre.euNOG_orthology.tab --peptides proteomes/chlre/uniprot-proteome%3AUP000006906_peptides.csv --output_dir proteomes/chlre/

#python scripts/define_grouping.py --spec phatc --grouping_type euNOG --grouping eggnog_output/phatc.euNOG_orthology.tab --peptides proteomes/phatc/phatc_uniprot_proteome_peptides.csv --output_dir proteomes/phatc/

#python scripts/define_grouping.py --spec lytva --grouping_type euNOG --grouping eggnog_output/lytva.euNOG_orthology.tab --peptides proteomes/lytva/Lvar2_2_peptides.csv --output_dir proteomes/lytva/
#python scripts/define_grouping.py --spec chqui --grouping_type protein  --peptides proteomes/chqui/Cquinoa_392_v1.0_peptides.csv --output_dir proteomes/chqui/

#python scripts/define_grouping.py --spec tetts --grouping_type protein  --peptides proteomes/tetts/uniprot_proteome_tetts_UP000009168_peptides.csv --output_dir proteomes/tetts/

#python scripts/define_grouping.py --spec tetts --grouping_type euNOG --grouping eggnog_mapper/uniprot_proteome_tetts_UP000009168.fasta.emapper.annotations.euNOGmapping --peptides proteomes/tetts/uniprot_proteome_tetts_UP000009168_peptides.csv --output_dir proteomes/tetts/

#python scripts/define_grouping.py --spec sollc --grouping_type protein  --peptides proteomes/sollc/sollc_uniprot_peptides.csv --output_dir proteomes/sollc/

#python scripts/define_grouping.py --spec sollc --grouping_type euNOG --grouping eggnog_mapper/sollc_uniprot.fasta.emapper.annotations.euNOGmapping  --peptides proteomes/sollc/sollc_uniprot_peptides.csv --output_dir proteomes/sollc/


#python scripts/define_grouping.py --spec cerri --grouping_type euNOG --grouping eggnog_output/cerri.euNOG_orthology.tab --peptides proteomes/cerri/McWhite201602_combined_prot_full_NR_uc097seeds_peptides.csv --output_dir proteomes/cerri/

#nohup python scripts/remove_nonidentifying_peptides.py --spec arath --grouping_type  protein --peptides proteomes/arath/uniprot-proteome%3AUP000006548_peptides.csv --output_dir proteomes/arath/ &
#nohup python scripts/remove_nonidentifying_peptides.py --spec orysj --grouping_type protein --peptides proteomes/orysj/uniprot-proteome%3AUP000000763_peptides.csv --output_dir proteomes/orysj/ &
#nohup python scripts/remove_nonidentifying_peptides.py --spec traes --grouping_type  protein --peptides proteomes/traes/uniprot-triticum%3AUP000019116_peptides.csv --output_dir proteomes/traes/ &
#nohup python scripts/remove_nonidentifying_peptides.py --spec braol --grouping_type  protein   --peptides proteomes/braol/uniprot-brassica+oleracea_peptides.csv --output_dir proteomes/braol/ &
#nohup python scripts/remove_nonidentifying_peptides.py --spec selml --grouping_type  protein  --peptides proteomes/selml/uniprot-proteome%3AUP000001514_peptides.csv --output_dir proteomes/selml/ &
#nohup python scripts/remove_nonidentifying_peptides.py --spec chlre --grouping_type  protein  --peptides proteomes/chlre/uniprot-proteome%3AUP000006906_peptides.csv --output_dir proteomes/chlre/ &
#nohup python scripts/remove_nonidentifying_peptides.py --spec cerri --grouping_type  protein --peptides proteomes/cerri/McWhite201602_combined_prot_full_NR_uc097seeds_peptides.csv --output_dir proteomes/cerri/ &

#python scripts/remove_nonidentifying_peptides.py --spec euggr --grouping_type protein  --peptides proteomes/euggr/Euglena_nr_peptides.csv --output_dir proteomes/euggr/
#python scripts/remove_nonidentifying_peptides2.py --spec phatc   --peptides proteomes/phatc/phatc_uniprot_proteome_peptides.csv --output_dir proteomes/phatc/
#remove_nonidentifying_peptides2.py

#python scripts/remove_nonidentifying_peptides2.py --spec lytva   --peptides proteomes/lytva/Lvar2_2_peptides.csv --output_dir proteomes/lytva/


############Lookup proteins/groups according to custom grouping or no grouping

#
#
#mkdir identified_elutions/sollc
#
#for exp in OP_Tomato_leaf_IEX2_201270607
#do
#    echo protein lookup $exp
#
#   python scripts/get_elution_profiles.py sollc protein $exp elutions/${exp}_elution.csv proteomes/sollc/unique_peptides_sollc_protein.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_prot.py identified_elutions/sollc/${exp}_elution_sollc_protein.csv
#
#   echo grouped lookup $exp
#   python scripts/get_elution_profiles.py sollc euNOG $exp elutions/${exp}_elution.csv proteomes/sollc/unique_peptides_sollc_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#
#   python scripts/get_wideform_group.py identified_elutions/sollc/${exp}_elution_sollc_euNOG.csv eggnog_mapper/sollc_uniprot.fasta.emapper.annotations.euNOGmapping annotation_files/all_annotations.csv
#
#
#done
#
#
mkdir identified_elutions/human
#

for exp in Anna_HEK293T_IEX_control_011217 Anna_HEK293T_IEX_mixed_bed_control_042217 Anna_HEK293T_IEX_mixed_bed_RNASE_050117 Anna_HEK293T_IEX_RNaseA_012517 Anna_HEK293T_SEC_Benzonase_121316 Anna_HEK293T_SEC_control_113016  Anna_HEK293T_SEC_RNaseA_120916 Anna_HEK293T_IEX_RNaseA_012517 RBC_1086_WWC_ghosts_prep1_11232016
do
echo grouped lookup $exp
   python scripts/get_elution_profiles.py human euNOG $exp elutions/${exp}_elution.csv proteomes/human/unique_peptides_human_euNOG.csv proteomes/contam/contam_benzo_peptides.csv

   python scripts/get_wideform_group.py identified_elutions/human/${exp}_elution_human_euNOG.csv eggnog_mapper/human_hmmer_euk.mapping annotation_files/all_annotations.csv


done



for exp in Hs_helaC_ph_tcs269_2 Hs_helaC_ph_tcs269_1 Hs_helaC_mar_SGF Hs_hekN_1201 Hs_helaN_ph_tcs375_P2 Hs_hekN_ph_hcw120_1 Hs_helaN_ph_tcs375_P3 Hs_helaN_mar_ief3_10 Hs_helaN_mar_ief5_8 Hs_CB660_1105 Hs_helaN_mar_SGF Hs_G166_1104 Hs_helaN_ph_hcw120_1 Hs_G166_1105 Hs_helaN_ph_hcw120_2 Hs_hek_kinase95 Hs_helaN_ph_saf48 Hs_hekN_1108 Hs_helaN_ph_tcs375_P1 Hs_hekN_ph_hcw120_2 Hs_helaN_ph_tcs375_P4 Hs_hek_wan_membraneopioid72 Hs_helaN_wan_hcw120_3 Hs_helaC_mar_ief3_10 Hs_wan_aguhplc108_1204 Hs_helaC_mar_ief3_10_repeat Hs_wan_shuhplc94_1206 Hs_helaC_mar_ief5_8 Hs_helaC_ph_wax43 Hs_helaN_1003 Hs_helaN_1010 

do
   echo protein lookup $exp

   python scripts/get_elution_profiles.py human protein $exp elutions/${exp}_elution.csv proteomes/human/unique_peptides_human_protein.csv proteomes/contam/contam_benzo_peptides.csv
   python scripts/get_wideform_prot.py identified_elutions/human/${exp}_elution_human_protein.csv

    echo group lookup $exp

   python scripts/get_elution_profiles.py human euNOG $exp elutions/${exp}_elution.csv proteomes/human/unique_peptides_human_euNOG.csv proteomes/contam/contam_benzo_peptides.csv

   python scripts/get_wideform_group.py identified_elutions/human/${exp}_elution_human_euNOG.csv eggnog_mapper/human_hmmer_euk.mapping annotation_files/all_annotations.csv




done



for exp in Ce_1104 Ce_6mg_1202 Ce_6mg_1203 Ce_beadsA_1206 Ce_beadsB_wan_60_1208 Ce_beadsC_wan_60_1209 Ce_beadsflow_1206 Ce_beadsL_1206 Ce_beadsN_wan_60_1208 Ce_beadsR_wan_60_1209 Ce_BNF_wan_60_1209 Ce_CRF_wan_60_1209 Ce_ups_1205 Ce_ups_1206 
do 

#   echo protein lookup $exp

#   python scripts/get_elution_profiles.py caeel protein $exp elutions/${exp}_elution.csv proteomes/caeel/unique_peptides_caeel_protein.csv proteomes/contam/contam_benzo_peptides.csv
   #python scripts/get_wideform_prot.py identified_elutions/caeel/${exp}_elution_caeel_protein.csv

   echo group lookup $exp

   python scripts/get_elution_profiles.py caeel euNOG $exp elutions/${exp}_elution.csv proteomes/caeel/unique_peptides_caeel_euNOG.csv proteomes/contam/contam_benzo_peptides.csv

   python scripts/get_wideform_group.py identified_elutions/caeel/${exp}_elution_caeel_euNOG.csv eggnog_mapper/caeel_hmmer_euk.mapping annotation_files/all_annotations.csv




done






for exp in Mm_1106 Mm_1201
do 

   echo protein lookup $exp

   python scripts/get_elution_profiles.py mouse protein $exp elutions/${exp}_elution.csv proteomes/mouse/unique_peptides_mouse_protein.csv proteomes/contam/contam_benzo_peptides.csv
   python scripts/get_wideform_prot.py identified_elutions/mouse/${exp}_elution_mouse_protein.csv


   echo group lookup $exp

   python scripts/get_elution_profiles.py mouse euNOG $exp elutions/${exp}_elution.csv proteomes/mouse/unique_peptides_mouse_euNOG.csv proteomes/contam/contam_benzo_peptides.csv

   python scripts/get_wideform_group.py identified_elutions/mouse/${exp}_elution_mouse_euNOG.csv eggnog_mapper/mouse_hmmer_euk.mapping annotation_files/all_annotations.csv




done

for exp in Dd_1110 Dd_1112
do 

   echo protein lookup $exp

   python scripts/get_elution_profiles.py dicdi protein $exp elutions/${exp}_elution.csv proteomes/dicdi/unique_peptides_dicdi_protein.csv proteomes/contam/contam_benzo_peptides.csv
   python scripts/get_wideform_prot.py identified_elutions/dicdi/${exp}_elution_dicdi_protein.csv

   echo group lookup $exp

   python scripts/get_elution_profiles.py dicdi euNOG $exp elutions/${exp}_elution.csv proteomes/dicdi/unique_peptides_dicdi_euNOG.csv proteomes/contam/contam_benzo_peptides.csv

   python scripts/get_wideform_group.py identified_elutions/dicdi/${exp}_elution_dicdi_euNOG.csv eggnog_mapper/dicdi_hmmer_euk.mapping annotation_files/all_annotations.csv




done

for exp in Dm_1107_Dem Dm_1107_DSL2 Dm_1201_Dem 
do
   echo protein lookup $exp

   python scripts/get_elution_profiles.py drome protein $exp elutions/${exp}_elution.csv proteomes/drome/unique_peptides_drome_protein.csv proteomes/contam/contam_benzo_peptides.csv
   python scripts/get_wideform_prot.py identified_elutions/drome/${exp}_elution_drome_protein.csv


   echo group lookup $exp

   python scripts/get_elution_profiles.py drome euNOG $exp elutions/${exp}_elution.csv proteomes/drome/unique_peptides_drome_euNOG.csv proteomes/contam/contam_benzo_peptides.csv

   python scripts/get_wideform_group.py identified_elutions/drome/${exp}_elution_drome_euNOG.csv eggnog_mapper/drome_hmmer_euk.mapping annotation_files/all_annotations.csv

done


for exp in Xl_muscrna_op_24_1302 Xl_caps_1308 Xl_musc_op_26_1212 Xl_musc_op_22_1303 
do 


   echo protein lookup $exp

   python scripts/get_elution_profiles.py xenla protein $exp elutions/${exp}_elution.csv proteomes/xenla/working_proteome/unique_peptides_xenla_protein.csv proteomes/contam/contam_benzo_peptides.csv
   python scripts/get_wideform_prot.py identified_elutions/xenla/${exp}_elution_xenla_protein.csv

   echo group lookup $exp

   python scripts/get_elution_profiles.py xenla euNOG $exp elutions/${exp}_elution.csv proteomes/xenla/working_proteome/unique_peptides_xenla_euNOG.csv proteomes/contam/contam_benzo_peptides.csv

   python scripts/get_wideform_group.py identified_elutions/xenla/${exp}_elution_xenla_euNOG.csv eggnog_mapper/xenla_hmmer_euk.mapping annotation_files/all_annotations.csv


done



for exp in OP_Tetrahymena_bodies_IEX_20170619 #OP_Tetrahymena_cilia_IEX_20170609
do

for exp in OP_Tetrahymena_cilia_SEC_20170708
do


    echo protein lookup $exp

   python scripts/get_elution_profiles.py tetts protein $exp elutions/${exp}_elution.csv proteomes/tetts/unique_peptides_tetts_protein.csv proteomes/contam/contam_benzo_peptides.csv
   python scripts/get_wideform_prot.py identified_elutions/tetts/${exp}_elution_tetts_protein.csv

   echo grouped lookup $exp
   python scripts/get_elution_profiles.py tetts euNOG $exp elutions/${exp}_elution.csv proteomes/tetts/unique_peptides_tetts_euNOG.csv proteomes/contam/contam_benzo_peptides.csv

   python scripts/get_wideform_group.py identified_elutions/tetts/${exp}_elution_tetts_euNOG.csv eggnog_mapper/tetts_diamond_euk.emapper.annotations.euNOG.mapping annotation_files/all_annotations.csv


done


#mkdir identified_elutions/lytva

#1 lytva
#for exp in OP_SeaUrchinSperm_WWC_20160119 
#do

   #protein lookup

#a   python scripts/get_elution_profiles.py lytva protein $exp elutions/${exp}_elution.csv proteomes/lytva/protein_unique_peptides_lytva.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_prot.py identified_elutions/lytva/${exp}_elution_lytva_protein.csv

#   grouped lookup
#   python scripts/get_elution_profiles.py lytva euNOG $exp elutions/${exp}_elution.csv proteomes/lytva/orthogroup_unique_peptides_lytva_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_group.py identified_elutions/lytva/${exp}_elution_lytva_euNOG.csv eggnog_output/lytva.euNOG_orthology.tab

#done

#mkdir identified_elutions/chqui

#1 chqui
#for exp in OP_Quinoa_IEX_20170513
#do
#
#   echo 'protein lookup'
#   python scripts/get_elution_profiles.py chqui protein $exp elutions/${exp}_elution.csv proteomes/chqui/unique_peptides_chqui_protein.csv proteomes/contam/contam_benzo_peptides.csv
#
#   echo 'Get wideform'
#   python scripts/get_wideform_prot.py identified_elutions/chqui/${exp}_elution_chqui_protein.csv
#
#   #grouped lookup
#   #python scripts/get_elution_profiles.py chqui euNOG $exp elutions/${exp}_elution.csv proteomes/chqui/orthogroup_unique_peptides_chqui_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#   #python scripts/get_wideform_group.py identified_elutions/chqui/${exp}_elution_chqui_euNOG.csv eggnog_output/chqui.euNOG_orthology.tab
#
#done



#mkdir identified_elutions/phatc

#2 phatc
#for exp in OP_Ptricornutum_SEC_20170403 OP_Ptricornutum_IEX_20170328  
#do

   #protein lookup

#   python scripts/get_elution_profiles.py phatc protein $exp elutions/${exp}_elution.csv proteomes/phatc/protein_unique_peptides_phatc.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_prot.py identified_elutions/phatc/${exp}_elution_phatc_protein.csv
#
#   grouped lookup
#   python scripts/get_elution_profiles.py phatc euNOG $exp elutions/${exp}_elution.csv proteomes/phatc/orthogroup_unique_peptides_phatc_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_group.py identified_elutions/phatc/${exp}_elution_phatc_euNOG.csv eggnog_output/phatc.euNOG_orthology.tab

#done





#mkdir identified_elutions/euggr

#1 euglena
#for exp in euglena_sec_1-2017  
#do

   #protein lookup

#   python scripts/get_elution_profiles.py euggr protein $exp elutions/${exp}_elution.csv proteomes/euggr/protein_unique_peptides_euggr.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_prot.py identified_elutions/euggr/${exp}_elution_euggr_protein.csv

   #grouped lookup
#   python scripts/get_elution_profiles.py euggr eunog $exp elutions/${exp}_elution.csv proteomes/euggr/orthogroup_unique_peptides_euggr_eunog.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_group.py identified_elutions/euggr/${exp}_elution_euggr_eunog.csv eggnog_output/euggr.euglena_nr.eunog_orthology.tab

#done

#exit 0

#5 Arabidopsis
#for exp in OP_Arabidopsis_sproutsSEC_02102017 At_Col_0_indark_201505 At_Col_0_indark_fraction_201504 At_Col_0_leaf_fraction_2014 At_Col_0_leaf_fraction_2015

#for exp in Arabidopsis_sprouts_WWC_2017 
#
#for exp in OP_Arabidopsis_mutantseed_20170323 
#do
#   #Protein lookup
#   python scripts/get_elution_profiles.py arath protein $exp elutions/${exp}_elution.csv proteomes/arath/protein_unique_peptides_arath.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_prot.py identified_elutions/arath/${exp}_elution_arath_protein.csv
#   #Grouped lookup
#   python scripts/get_elution_profiles.py arath euNOG $exp elutions/${exp}_elution.csv proteomes/arath/orthogroup_unique_peptides_arath_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_group.py identified_elutions/arath/${exp}_elution_arath_euNOG.csv eggnog_output/arath.euNOG_orthology.tab
#done
#
#
#
#
#
###3 Rice
#for exp in RiceLeaf_SEC_04292016 Rice_201505_uniprot RiceL_IEX
#do
#   #Protein lookup
#   python scripts/get_elution_profiles.py orysj protein $exp elutions/${exp}_elution.csv proteomes/orysj/protein_unique_peptides_orysj.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_prot.py identified_elutions/orysj/${exp}_elution_orysj_protein.csv
#   #Grouped lookup
#   python scripts/get_elution_profiles.py orysj euNOG $exp elutions/${exp}_elution.csv proteomes/orysj/orthogroup_unique_peptides_orysj_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_group.py identified_elutions/orysj/${exp}_elution_orysj_euNOG.csv eggnog_output/orysj.euNOG_orthology.tab
#done
#
##5 Wheat
#for exp in OP_WheatGrassSEC_20170202 wgIEF wgSEC1 wheatgermIEX WheatGermSEC_07-2015 
#do
#   #Protein lookup
#   python scripts/get_elution_profiles.py traes protein $exp elutions/${exp}_elution.csv proteomes/traes/protein_unique_peptides_traes.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_prot.py identified_elutions/traes/${exp}_elution_traes_protein.csv
#   #Grouped lookup
#   python scripts/get_elution_profiles.py traes euNOG $exp elutions/${exp}_elution.csv proteomes/traes/orthogroup_unique_peptides_traes_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#   python scripts/get_wideform_group.py identified_elutions/traes/${exp}_elution_traes_euNOG.csv eggnog_output/traes.euNOG_orthology.tab
#done
#
##3 Broccoli
#
#
#for exp in BroccoliNE_WWC Broccolinuclei_6-2016 Broccolinuclei_IEF_9-2016
#for exp in OP_Broccoli_leaf_IEX2_20170317 
#do
#    #Protein lookup
#    python scripts/get_elution_profiles.py braol protein $exp elutions/${exp}_elution.csv proteomes/braol/protein_unique_peptides_braol.csv proteomes/contam/contam_benzo_peptides.csv
#    python scripts/get_wideform_prot.py identified_elutions/braol/${exp}_elution_braol_protein.csv
#    #Grouped lookup
#    python scripts/get_elution_profiles.py braol euNOG $exp elutions/${exp}_elution.csv proteomes/braol/orthogroup_unique_peptides_braol_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#    python scripts/get_wideform_group.py identified_elutions/braol/${exp}_elution_braol_euNOG.csv eggnog_output/braol.euNOG_orthology.tab
#done
#
##2 Selaginella
#for exp in OP_SelaginellaSEC_20160309 selaginella_WWC
#do
#    #Protein lookup
#    python scripts/get_elution_profiles.py selml protein $exp elutions/${exp}_elution.csv proteomes/selml/protein_unique_peptides_selml.csv proteomes/contam/contam_benzo_peptides.csv
#    python scripts/get_wideform_prot.py identified_elutions/selml/${exp}_elution_selml_protein.csv
#    #Grouped lookup
#    python scripts/get_elution_profiles.py selml euNOG $exp elutions/${exp}_elution.csv proteomes/selml/orthogroup_unique_peptides_selml_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#    python scripts/get_wideform_group.py identified_elutions/selml/${exp}_elution_selml_euNOG.csv eggnog_output/selml.euNOG_orthology.tab
#done
#
##2 Chlamydomonas
#for exp in OP_Chlamydomonas_SEC Chlamydomonas_WWC_9-2016 OP_Chlamy_Nuclei_SEC_20170428
#for exp in OP_Chlamy_Nuclei_SEC_20170428
#do
#    #Protein lookup
#    python scripts/get_elution_profiles.py chlre protein $exp elutions/${exp}_elution.csv proteomes/chlre/protein_unique_peptides_chlre.csv proteomes/contam/contam_benzo_peptides.csv
#    python scripts/get_wideform_prot.py identified_elutions/chlre/${exp}_elution_chlre_protein.csv
#    #Grouped lookup
#    python scripts/get_elution_profiles.py chlre euNOG $exp elutions/${exp}_elution.csv proteomes/chlre/orthogroup_unique_peptides_chlre_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#    python scripts/get_wideform_group.py identified_elutions/chlre/${exp}_elution_chlre_euNOG.csv eggnog_output/chlre.euNOG_orthology.tab
#done
#
#
##2 Fern
#for exp in Fern_Frond_WWC_20151217-20160119 OP_Fernfrond_SEC_112016
#do
#    #Protein lookup
#    python scripts/get_elution_profiles.py cerri protein $exp elutions/${exp}_elution.csv proteomes/cerri/protein_unique_peptides_cerri.csv proteomes/contam/contam_benzo_peptides.csv
#    python scripts/get_wideform_prot.py identified_elutions/cerri/${exp}_elution_cerri_protein.csv
#    #Grouped lookup
#    python scripts/get_elution_profiles.py cerri euNOG $exp elutions/${exp}_elution.csv proteomes/cerri/orthogroup_unique_peptides_cerri_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
#    python scripts/get_wideform_group.py identified_elutions/cerri/${exp}_elution_cerri_euNOG.csv eggnog_output/cerri.euNOG_orthology.tab
#done
#
#
#
#
#
#
##for exp in Hs_CB660_1105 Hs_hekN_ph_hcw120_2 Hs_helaC_ph_tcs269_2 Hs_helaN_ph_hcw120_1 Hs_helaN_wan_hcw120_3 Hs_G166_1104 Hs_hek_wan_membraneopioid72 Hs_helaC_ph_wax43  Hs_helaN_ph_hcw120_2 Hs_wan_aguhplc108_1204 Hs_G166_1105 Hs_helaC_mar_ief3_10 Hs_helaN_1003  Hs_helaN_ph_saf48  Hs_wan_shuhplc94_1206 Hs_hek_kinase95 Hs_helaC_mar_ief3_10_repeat Hs_helaN_1010  Hs_helaN_ph_tcs375_P1 Hs_hekN_1108 Hs_helaC_mar_ief5_8 Hs_helaN_mar_ief3_10 Hs_helaN_ph_tcs375_P2 Hs_hekN_1201 Hs_helaC_mar_SGF Hs_helaN_mar_ief5_8 Hs_helaN_ph_tcs375_P3 Hs_hekN_ph_hcw120_1 Hs_helaC_ph_tcs269_1 Hs_helaN_mar_SGF  Hs_helaN_ph_tcs375_P4 
##do
##
##    python scripts/get_elution_profiles.py human protein $exp elutions/${exp}_elution.csv proteomes/human/uniprot-proteome_human_reviewed_peptides.csv proteomes/contam/contam_benzo_peptides.csv
##    python scripts/get_wideform_prot.py identified_elutions/human/${exp}_elution_human_protein.csv
##done
##
##
##python scripts/get_elution_profiles.py human euNOG Hs_hekN_1108_elution elutions/Hs_hekN_1108 proteomes/human/orthogroup_unique_peptides_human_euNOG.csv proteomes/contam/contam_benzo_peptides.csv
##
#python scripts/get_elution_profiles.py human protein Hs_hekN_1108_elution elutions/Hs_hekN_1108_elution.csv proteomes/human/uniprot-proteome_human_reviewed_peptides.csv proteomes/contam/contam_benzo_peptides.csv
#
#python scripts/get_wideform_prot.py identified_elutions/human/Hs_hekN_1108_elution_elution_human_protein.csv
#
