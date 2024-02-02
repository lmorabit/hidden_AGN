rule generateZmaxes:
	input:
		"src/static/lockman_final_cross_match_catalogue-v1.0.fits",
		"src/static/lockman_rms_starmask_optical.fits"
	output:
		"src/data/RLF.fits",
		"src/data/zmaxes.fits"
	conda:
		"environment.yml"
	script:
		"src/scripts/generateZmaxes.py"

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
                "src/data/RLF.fits"
        output:
                "mauch_sadler_RLFs.png"
        conda:
                "environment.yml"
        script:
                "src/scripts/mauchsadler.py"

