# TE-Networks
Modification of scripts from the https://github.com/MiguelMSandin/SSNetworks that works better for fragmented hits, such as Transposable Elements


The modified scripts are:

**1.1_blastn_allAgainstAll.sh**, which does an all-vs-all blast of the sequences in a fasta file. It is modified from Miguel's scrips to use blastn instead of megablast


**1.2_blastClean_modIW_v2.py**, which cleans up the blast output by removing self-hits and reciprocal hits. In my modified version it also deals with fragmented and overlapping hits which you get with blastn of transposable elements. 

**2.1_buildNetwork_IWmod**, which builds the network by considering coverage and percent identity. It is modified to calculate it from the output of the modified scrips. 
