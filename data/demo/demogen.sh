#!/bin/bash
# header
head -1 /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep1_sites.bed > FrontalCortex_rep1_sites.bed
head -1 /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep1_sites.bed > FrontalCortex_rep2_sites.bed
head -1 /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep1_sites.bed > Liver_rep1_sites.bed
head -1 /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep1_sites.bed > Liver_rep2_sites.bed
head -1 /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep1_sites.bed > HeartF_rep1_sites.bed
head -1 /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep1_sites.bed > HeartF_rep2_sites.bed
# ZNF609
grep "64791491\|64792365" /data/circrna/Human/ENCODE_2014/FrontalCortex/FrontalCortex_rep1_sites.bed | grep -v 64915251 >> FrontalCortex_rep1_sites.bed
grep "64791491\|64792365" /data/circrna/Human/ENCODE_2014/FrontalCortex/FrontalCortex_rep2_sites.bed | grep -v 64915251 >> FrontalCortex_rep2_sites.bed

grep "64791491\|64792365" /data/circrna/Human/ENCODE_2014/Liver/Liver_rep1_sites.bed >> Liver_rep1_sites.bed
grep "64791491\|64792365" /data/circrna/Human/ENCODE_2014/Liver/Liver_rep2_sites.bed | grep -v 64802316 >> Liver_rep2_sites.bed

grep "64791491\|64792365" /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep1_sites.bed >> HeartF_rep1_sites.bed
grep "64791491\|64792365" /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep2_sites.bed >> HeartF_rep2_sites.bed

# SLC8A1

grep "40655612\|40657441" /data/circrna/Human/ENCODE_2014/FrontalCortex/FrontalCortex_rep1_sites.bed | grep -v "40659437\|40505758\|40673788" >> FrontalCortex_rep1_sites.bed
grep "40655612\|40657441" /data/circrna/Human/ENCODE_2014/FrontalCortex/FrontalCortex_rep2_sites.bed | grep -v 40657110 >> FrontalCortex_rep2_sites.bed

grep "40655612\|40657441" /data/circrna/Human/ENCODE_2014/Liver/Liver_rep1_sites.bed >> Liver_rep1_sites.bed
grep "40655612\|40657441" /data/circrna/Human/ENCODE_2014/Liver/Liver_rep2_sites.bed | grep -v 64802316 >> Liver_rep2_sites.bed

grep "40655612\|40657441" /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep1_sites.bed | grep -v "40602403\|40659437\|40511939\|40673788\|40549690\|40655608\|40709870\|40710964\|40657308\|40696663\|40729232\|40680489\|40604814\|40476069\|40698853\|40729228\|40653214" >> HeartF_rep1_sites.bed
grep "40655612\|40657441" /data/circrna/Human/ENCODE_2014/Heart_fetal/HeartF_rep2_sites.bed | grep -v "40602403\|40659437\|40511939\|40673788\|40549690\|40655608\|40709870\|40710964\|40657308\|40696663\|40729232\|40680489\|40604814\|40476069\|40698853\|40729228\|40653214" >> HeartF_rep2_sites.bed
