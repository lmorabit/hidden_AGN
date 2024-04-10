rule generateVmaxes:
	input:
		"src/static/lockman_final_cross_match_catalogue-v1.0_classifications_catalogue_filtered_full_SNR5_fluxscaled_withoffset_noduplicates_with_lotss_DR1_detectable.fits",
		"src/static/lockman_rms_starmask_optical.fits",
                "src/static/cochrane_2023_tableA1.csv",
                "src/static/kondapally_2022_table1.csv"
	output:
		"src/data/lockman_6arcsec_vmax.fits"
	conda:
		"environment.yml"
	script:
		"src/scripts/generateVmaxes.py"

rule testplot:
        input:
                "src/static/RLFS_50MYR_SF.csv",
                "src/static/RLFS_50MYR_AGN.csv"
        output:
                "test.png"
        conda:
                "environment.yml"
        script:
                "src/scripts/testplot.py"
rule completeness:
	input:
		"src/static/cochrane_2023_table1.csv",
		"src/static/cochrane_2023_tableA1.csv",
		"src/static/kondapally_2022_table1.csv",
		"src/static/kondapally_2022_table2.csv"
	output:
		"completeness.png"
	conda:
		"environment.yml"
	script:
		"src/scripts/completeness.py"
rule mauchsadler:
        input:
                "src/static/mauch_sadler_table5.csv",
                "src/data/lockman_vmaxes.fits",
		"src/static/kondapally_2022_table2.csv",
		"src/static/cochrane_2023_table1.csv"
        output:
                "mauch_sadler_RLFs.png"
        conda:
                "environment.yml"
        script:
                "src/scripts/mauchsadler.py"

