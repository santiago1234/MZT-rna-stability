import pandas as pd
import os
from optimalcodon.limsdatacopy import get_fastqs_from_lims

path_to_fastqs = "/n/analysis/Bazzini/arb/MOLNG-2541/HGYJ3BGX9/" # TODO: change this to unix
destination_dir = "rawfastqfiles"
sample_report = os.path.join(path_to_fastqs, "Sample_Report.csv")
sample_report = pd.read_csv(sample_report)

# subset the data for the given experiment
MOLNG_2541 = sample_report[(sample_report.Reference == "danRer10") & (sample_report.Order == "MOLNG-2541")]
MOLNG_2541 = MOLNG_2541.drop_duplicates(['IndexSequence1', 'SampleName'])

for _, row in MOLNG_2541.iterrows():
    destinationname = os.path.join(destination_dir, 'zfishRibo0' + row.SampleName + ".fastq.gz")
    get_fastqs_from_lims(path_to_fastqs, row.IndexSequence1, destinationname)




### MOLNG-2540
path_to_fastqs = "/n/analysis/Bazzini/arb/MOLNG-2540/HFVYFBGX9/" # TODO: change this to unix
destination_dir = "rawfastqfiles"
sample_report = os.path.join(path_to_fastqs, "Sample_Report.csv")
sample_report = pd.read_csv(sample_report)

# subset the data for the given experiment
MOLNG_2540 = sample_report[(sample_report.Reference == "danRer10") & (sample_report.Order == "MOLNG-2540")]
MOLNG_2540 = MOLNG_2540.drop_duplicates(['IndexSequence1', 'SampleName'])

for _, row in MOLNG_2540.iterrows():
    destinationname = os.path.join(destination_dir, 'zfishPolyA' + row.SampleName + ".fastq.gz")
    get_fastqs_from_lims(path_to_fastqs, row.IndexSequence1, destinationname)


# MOLNG-2539: Zebrafish_alpha_timecourse
path_to_fastqs = "/n/analysis/Bazzini/arb/MOLNG-2539/HFWJMBGX9/" # TODO: change this to unix
destination_dir = "rawfastqfiles"
sample_report = os.path.join(path_to_fastqs, "Sample_Report.csv")
sample_report = pd.read_csv(sample_report)

# subset the data for the given experiment
MOLNG_2539 = sample_report[(sample_report.Reference == "danRer10") & (sample_report.Order == "MOLNG-2539")]
MOLNG_2539 = MOLNG_2539.drop_duplicates(['IndexSequence1', 'SampleName'])

for _, row in MOLNG_2539.iterrows():
    destinationname = os.path.join(destination_dir, 'TreatedAamanitin_zfishPolyA' + row.SampleName + ".fastq.gz")
    get_fastqs_from_lims(path_to_fastqs, row.IndexSequence1, destinationname)

