#!/bin/bash
#hostname

grep_list=("1B4B" "1B4E" "1H7N" "1H7R" "1II8" "1L6S" "1PV8" "1SXJ" "1W1W" "1W4R" "1W5Q" "1XX6" "1XXA" "2B8T" "2C14" "2O71" "2OF5" "2P5M" "2P5T" "2ZCI" "2ZFZ" "3AUY" "3CAG" "3CZP" "3DU6" "3EZQ" "3O8O" "3OBK" "3OQ9" "3QKS" "4AG6" "4EW5" "4HX0" "4OL8" "4UXH" "4XGC" "5CQG" "5FMF" "5G2X" "5HHJ" "5HHL" "5HMS" "5JVO" "5LZL" "5UJ7" "5Z69" "6AC0" "6AR1" "6AR3" "6ME0" "6MSB")
for pdbid in $grep_list; do
	grep $pdbid ./bc-70.out >>greped_pdbids.txt
done

TN_list=("1A1D" "1AAR" "1AW5" "1B4B" "1B4E" "1B4K" "1C3T" "1C9B" "1CDW" "1CMX" "1D3Z" "1DDF" "1DZF" "1E3Y" "1E41" "1E51" "1EB3" "1F9J" "1FXT" "1G6J" "1GJP" "1GZG" "1H7N" "1H7O" "1H7P" "1H7R" "1I3Q" "1I50" "1I6H" "1I8J" "1II8" "1IRU" "1IYJ" "1JFI" "1K83" "1L6S" "1L6Y" "1MIU" "1MJE" "1NBF" "1NGM" "1NH2" "1NIK" "1NT9" "1NVP" "1OGW" "1OHL" "1OTR" "1P3Q" "1PLQ" "1PLR" "1PQV" "1PV8" "1Q0W" "1Q5W" "1QML" "1QN3" "1QN4" "1QN5" "1QN6" "1QN7" "1QN8" "1QN9" "1QNA" "1QNB" "1QNC" "1QNE" "1QNV" "1R5U" "1R9S" "1R9T" "1RM1" "1S1Q" "1SFO" "1SIF" "1SXJ" "1TBA" "1TBE" "1TBP" "1TGH" "1TWA" "1TWC" "1TWF" "1TWG" "1TWH" "1UBI" "1UBQ" "1UD7" "1UZX" "1V80" "1V81" "1VOK" "1VOL" "1W1W" "1W31" "1W4R" "1W54" "1W56" "1W5M" "1W5N" "1W5O" "1W5P" "1W5Q" "1WCM" "1WR1" "1WR6" "1WRD" "1XBT" "1XD3" "1XQQ" "1XX6" "1XXA" "1XXB" "1XXC" "1Y14" "1Y1V" "1Y1W" "1Y1Y" "1Y77" "1YD8" "1YIW" "1YJ1" "1YLV" "1YTB" "1YTF" "1ZGU" "1ZW7" "2AYO" "2B63" "2B8K" "2B8T" "2BGF" "2C13" "2C14" "2C15" "2C16" "2C18" "2C19" "2C7M" "2C7N" "2D3G" "2DEN" "2DX5" "2E2H" "2E2I" "2E2J" "2EWR" "2FCL" "2FCM" "2FCN" "2FCQ" "2FCS" "2FID" "2FIF" "2FUH" "2G3Q" "2G45" "2GBJ" "2GBK" "2GBM" "2GBN" "2GBR" "2GMI" "2HD5" "2HTH" "2IBI" "2J7Q" "2JA5" "2JA6" "2JA7" "2JA8" "2JF5" "2JRI" "2JT4" "2JVC" "2JWZ" "2JY6" "2JZZ" "2K39" "2K6D" "2K8B" "2K8C" "2KDE" "2KDF" "2KHW" "2KJH" "2KLG" "2KN5" "2KOX" "2KTF" "2KWU" "2KWV" "2L00" "2L0F" "2L0T" "2L3Z" "2LD9" "2LJ5" "2LVO" "2LVP" "2LVQ" "2LZ6" "2M0X" "2MBB" "2MBH" "2MBO" "2MBQ" "2MCN" "2MI8" "2MJ5" "2MJB" "2MOR" "2MRE" "2MRO" "2MSG" "2MUR" "2MWS" "2N13" "2N2K" "2N3U" "2N3V" "2N3W" "2NBD" "2NBE" "2NR2" "2NVQ" "2NVT" "2NVX" "2NVY" "2NVZ" "2O6V" "2O71" "2OD8" "2OF5" "2OOB" "2P5M" "2P5T" "2PE9" "2PEA" "2QHO" "2R7Z" "2R92" "2R93" "2RR9" "2RSU" "2RU6" "2UZ3" "2VUM" "2WDT" "2WOQ" "2WVJ" "2WWZ" "2WX0" "2WX1" "2XBB" "2XEW" "2XK5" "2YU9" "2Z0I" "2Z1B" "2Z59" "2ZCB" "2ZCC" "2ZCI" "2ZFZ" "2ZNV" "3A1Q" "3A33" "3A9J" "3A9K" "3ALB" "3AUL" "3AUX" "3AUY" "3AV0" "3BUE" "3BY4" "3C0R" "3CAG" "3CMM" "3CQZ" "3CZP" "3DU5" "3DU6" "3DVG" "3DVN" "3EEC" "3EFU" "3EHV" "3EIK" "3EZQ" "3F1W" "3FKI" "3GPM" "3GPN" "3GTG" "3GTJ" "3GTK" "3GTL" "3GTM" "3GTO" "3GTP" "3GTQ" "3H0G" "3H1U" "3H3V" "3H7P" "3H7S" "3HM3" "3HOU" "3HOV" "3HOW" "3HOX" "3HOY" "3HOZ" "3I3T" "3I4M" "3I4N" "3IFW" "3IHP" "3J0K" "3J1N" "3JCO" "3JCP" "3JSV" "3JVZ" "3JW0" "3K1F" "3K7A" "3K9P" "3KVF" "3KW5" "3KYL" "3LDZ" "3M3J" "3M3Y" "3M4O" "3MHS" "3MTN" "3N30" "3N32" "3N3K" "3NHE" "3NOB" "3NS8" "3O65" "3O8O" "3OBK" "3OC3" "3OCI" "3OFI" "3OJ3" "3OJ4" "3OLM" "3ONS" "3OQ9" "3PHD" "3PHW" "3PO2" "3PO3" "3PRM" "3PRP" "3PT2" "3PTF" "3QKR" "3QKS" "3QKU" "3QT1" "3RUL" "3RZD" "3RZO" "3S14" "3S15" "3S16" "3S17" "3S1M" "3S1N" "3S1Q" "3S1R" "3S2D" "3S2H" "3T5X" "3TBL" "3TMP" "3TXM" "3TXN" "3UGB" "3UNB" "3UNE" "3UNF" "3UNH" "3V60" "3V61" "3V62" "3V6C" "3V6E" "3VFK" "3VHT" "3VUW" "3VUX" "3VUY" "3WWQ" "3WXG" "3ZLZ" "3ZNH" "3ZNI" "4A3B" "4A3C" "4A3D" "4A3E" "4A3F" "4A3G" "4A3I" "4A3J" "4A3K" "4A3L" "4A3M" "4A93" "4AG5" "4AG6" "4AP4" "4AUQ" "4BBN" "4BBR" "4BBS" "4BOS" "4BOZ" "4BVU" "4BXX" "4BXZ" "4BY1" "4BY7" "4C2M" "4C3H" "4C3I" "4C3J" "4CR2" "4CR3" "4CR4" "4DDG" "4DDI" "4DHJ" "4DHZ" "4EW5" "4FJV" "4GSW" "4GU2" "4HCN" "4HJK" "4HK2" "4HX0" "4HXD" "4I6L" "4I6N" "4IG7" "4II2" "4II3" "4IUM" "4JIO" "4JQW" "4K1R" "4K7S" "4K7U" "4K7W" "4KSK" "4L60" "4L6P" "4LCD" "4LDT" "4LJO" "4LJP" "4M0W" "4MDK" "4MM3" "4MSM" "4MSQ" "4NNJ" "4NQK" "4NQL" "4OL8" "4P4H" "4PIG" "4PIH" "4PIJ" "4PQT" "4Q5E" "4Q5H" "4R3O" "4R62" "4R67" "4RF0" "4RF1" "4ROC" "4ROD" "4ROE" "4S1Z" "4S22" "4UEL" "4UF6" "4UN2" "4UXH" "4UXI" "4UXJ" "4V1M" "4V1N" "4V1O" "4V3K" "4V3L" "4WHV" "4WLR" "4WUR" "4WZP" "4WZS" "4X67" "4X6A" "4XGC" "4XKH" "4XKL" "4XOF" "4XOK" "4XOL" "4XYZ" "4Y1H" "4Y52" "4Y7N" "4YHR" "4YM7" "4Z9S" "4ZFR" "4ZFT" "4ZPZ" "4ZUX" "5A0Q" "5A5B" "5AF4" "5AF5" "5AF6" "5AIT" "5AIU" "5BNB" "5BZ0" "5C3E" "5C44" "5C4A" "5C4J" "5C4X" "5C7J" "5C7M" "5CAW" "5CQG" "5CRA" "5CVM" "5CVN" "5CVO" "5D0K" "5D0M" "5DFL" "5DK8" "5DNY" "5DSV" "5E6J" "5EDV" "5EMZ" "5EYA" "5F3W" "5FER" "5FJ8" "5FJ9" "5FJA" "5FLM" "5FMF" "5FYW" "5FZ5" "5G2X" "5G5L" "5GJQ" "5GJR" "5GO7" "5GO8" "5GOB" "5GOC" "5GOD" "5GOG" "5GOH" "5GOI" "5GOJ" "5GOK" "5GVI" "5H7S" "5HHJ" "5HHK" "5HHL" "5HMS" "5HNR" "5HPK" "5HPL" "5HPS" "5HPT" "5IBK" "5IFR" "5IP7" "5IP9" "5IRF" "5IRG" "5IY6" "5IY7" "5IY8" "5IY9" "5IYA" "5IYB" "5IYC" "5IYD" "5J26" "5J8P" "5JBV" "5JBY" "5JG6" "5JNE" "5JP3" "5JQS" "5JTJ" "5JTV" "5JVO" "5K9P" "5KGF" "5KHY" "5KNL" "5KYC" "5KYD" "5KYE" "5KYF" "5L4G" "5L4K" "5L6H" "5L6I" "5L6J" "5L8H" "5L8W" "5LE5" "5LEX" "5LEY" "5LEZ" "5LF0" "5LF1" "5LF3" "5LF4" "5LF6" "5LF7" "5LMX" "5LN1" "5LN3" "5LRV" "5LRW" "5LRX" "5LZL" "5M32" "5M3F" "5M3M" "5M5W" "5M5X" "5M5Y" "5M64" "5M93" "5MHB" "5MN9" "5MNJ" "5MP9" "5MPA" "5MPB" "5MPC" "5N2W" "5N38" "5N5Y" "5N5Z" "5N60" "5N61" "5N9G" "5NL4" "5NL5" "5NLF" "5NLI" "5NLJ" "5NMC" "5NVG" "5O44" "5O6S" "5O6T" "5O76" "5OA1" "5OHK" "5OHL" "5OHM" "5OHN" "5OHP" "5OHV" "5OIK" "5OQJ" "5OQM" "5OT2" "5OXH" "5OXI" "5SVA" "5T0C" "5T0G" "5T0H" "5T0I" "5T0J" "5T9D" "5TOF" "5TOG" "5TR4" "5TTE" "5TUT" "5TXK" "5U0S" "5U4P" "5U5Q" "5UDH" "5UJ7" "5UJL" "5UJM" "5UJN" "5ULF" "5ULH" "5ULK" "5V1Y" "5V1Z" "5V5G" "5V5H" "5V5I" "5V69" "5V6A" "5V7K" "5V7L" "5V7M" "5VBT" "5VEY" "5VF0" "5VFO" "5VFP" "5VFQ" "5VFR" "5VFS" "5VFT" "5VFU" "5VGZ" "5VHF" "5VHH" "5VHI" "5VHM" "5VHN" "5VHO" "5VHQ" "5VHR" "5VHS" "5VNZ" "5VO0" "5VVR" "5VVS" "5VZM" "5VZW" "5W46" "5W4U" "5W51" "5W5Y" "5W64" "5W65" "5W66" "5WFI" "5WVI" "5WVK" "5WVO" "5X3M" "5X3N" "5X3O" "5X4Z" "5X50" "5X51" "5XBO" "5XDP" "5XIS" "5XIT" "5XIU" "5XK4" "5XK5" "5XME" "5XOG" "5XON" "5XPK" "5XU8" "5XVE" "5YDK" "5YDR" "5YIJ" "5YIK" "5YMY" "5YT6" "5Z67" "5Z68" "5Z69" "5ZBU" "5ZD0" "5ZQ3" "5ZQ4" "5ZQ5" "5ZQ6" "5ZQ7" "6A5L" "6A5O" "6A5P" "6A5R" "6A5T" "6A5U" "6A6I" "6AC0" "6AQR" "6AR1" "6AR3" "6ASR" "6AVO" "6B7M" "6B7O" "6BLO" "6BLP" "6BM2" "6BM4" "6BQF" "6BVA" "6BYH" "6C16" "6CNB" "6CNC" "6CND" "6CNF" "6CP2" "6CPM" "6CRN" "6CX2" "6CX3" "6CX4" "6D0Q" "6D0R" "6D4P" "6D68" "6D6I" "6DC6" "6DGF" "6DJ9" "6DJW" "6DJX" "6DRD" "6E2B" "6E49" "6E53" "6E5B" "6EF3" "6EI1" "6EPC" "6EPD" "6EPE" "6EPF" "6EQI" "6EU0" "6EU1" "6EU2" "6EU3" "6EXV" "6F40" "6F41" "6F42" "6F44" "6FDK" "6FGE" "6FTX" "6FVT" "6FVU" "6FVV" "6FVW" "6FVX" "6FVY" "6FX4" "6FYH" "6GLC" "6GMH" "6GML" "6GYK" "6GYL" "6GYM" "6GZS" "6H4H" "6H67" "6H68" "6HEI" "6HEK" "6HKO" "6HLQ" "6HLR" "6HLS" "6HPR" "6I84" "6IF1" "6INQ" "6IR9" "6ISU" "6J2C" "6J2N" "6J2Q" "6J2X" "6J30" "6J4W" "6J4X" "6J4Y" "6J4Z" "6J50" "6J51" "6J99" "6JB6" "6JB7" "6JKY" "6JMA" "6JWI" "6K4I" "6KIU" "6KIV" "6KIW" "6KOW" "6KOX" "6ME0" "6MEC" "6ML1" "6MSB" "6MSD" "6MSE" "6MSG" "6MSH" "6MSJ" "6MSK" "6N13" "6NJ9" "6NJG" "6NJQ" "6NN6" "6NOG" "6NQA" "6NYA" "6NYO" "6O6C" "6O82" "6O83" "6O96" "6O9L" "6OA9" "6OAA" "6OB1" "6OI4" "6OQ1" "6OQ2" "6PGV" "6QF8" "6QK9" "6R70" "6REY" "6RGQ" "6RQH" "6RQL" "6RQT" "6RRD" "6RUI" "6RUO" "6RWE" "6RYA" "6S53" "6UH5" "6V5D")
Merged_list=($TN_list $grep_list)

# Fetch chain sequences from PDB
p1="https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList="
p2="&compressionType=uncompressed"
for TN in $TN_list; do
	url2get="$p1$TN$p2"
	wget $url2get
done
cat down* >TNseqs.faa

# (try to) Find UniProtID for TNpdbIDs
for pdbid in $Merged_list; do
	grep $pdbid ./uniprot_pdb.tsv >>greped_uniprot_pdb.tsv
done
sort -u greped_uniprot_pdb.tsv >sorted_greped_uniprot_pdb.tsv
rm greped_uniprot_pdb.tsv
awk '{print $1}' sorted_greped_uniprot_pdb.tsv | sort -u >UniProtIDs2get.txt

# (try to) Find Pfams for TNpdbIDs
for pdbid in $Merged_list; do
	grep $pdbid ./pdb_chain_pfam.tsv >>greped_pfam_pdb.tsv
done
sort -u greped_pfam_pdb.tsv >sorted_greped_pfam_pdb.tsv
rm greped_pfam_pdb.tsv
awk '{print $4}' sorted_greped_pfam_pdb.tsv | sort -u >PFamsIDs2get.txt
grep -f PFamsIDs2get.txt /media/uri/HDD1/uri/DBs/Pfam32_A_info.tsv -A1 >greped_pfams.faa

grep -f PFamsIDs2get.txt $(find "rmdup_Pfam-A.fasta") -A1 >greped_pfams.faa
cat greped_pfams.faa rmdup_linear_TNseqs.faa >>TNseqs.faa
mkdir MMseqs2
cd MMseqs2
mmseqs2 easy-linclust ../TNseqs.faa TNs_seqs.clu tmp --min-seq-id 0.99 -c 1.0 --cov-mode 1 --kmer-per-seq-scale 0.4 --threads 7
mv ./TNs_seqs.clu_rep_seq.fasta ../TNseqs_reps.faa
cd ../
rm MMseqs2 -r
