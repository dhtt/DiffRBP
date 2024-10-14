# utilities.py>

from os import path
import subprocess

def get_tissues(matching_files, PPICompare_options):
	matching_df = pd.read_csv(matching_files, sep=" ", header=None)
	if "-d" in PPICompare_options:
		if "-m" in PPICompare_options:
			matching_df.columns = ["count_file", "ppin_file", "ddin_file", "major_transcripts_file"]
		else:
			matching_df.columns = ["count_file", "ppin_file", "ddin_file"]
	else:
		matching_df.columns = ["count_file", "ppin_file"]
    
	matching_df['tissue'] = matching_df['count_file'].apply(lambda x: x.split("/")[-1].split("_")[0])
	return matching_df

def move_file(old_dir: str, new_dir: str, *args):
	if not path.exists(new_dir):
		subprocess(f"mkdir -p {new_dir}")
	for file in args:
		if path.exists(f"{old_dir}/{file}"):
			subprocess(f"cp {old_dir}/{file} {new_dir}")
