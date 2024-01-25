rule comparison_plot:
        input:
                "src/data/kondapally_2022_table2.csv",
                "src/data/cochrane_2023_table1.csv",
                "src/data/rlfs"
        output:
                "deep_fields_RLFs.png"
        conda:
                "environment.yml"
        script:
                "src/scripts/comparison_plot.py"
rule rlf_evolution:
	input:
		"src/data/rlfs",
		"src/data/vmaxes"
	output:
		"src/tex/figures/RLF_evolution.png",
                "src/tex/output/integrated_differences.txt",
		"src/tex/output/average_integrated_differences.txt"
	conda:
		"environment.yml"
	script:
		"src/scripts/RLF_evolution.py"
rule calculate_vars:
	input: 
		"src/data/en1_03_matched_inMOC_inHR.fits",
		"src/data/lockman_03_matched_inMOC_inHR.fits",
	output:
                "src/tex/output/en1_detectable.txt",
                "src/tex/output/flowchart_numbers.txt",
                "src/tex/output/lockman_detectable.txt"
	conda:
		"environment.yml"
	script:
		"src/scripts/flowchart_numbers.py"
