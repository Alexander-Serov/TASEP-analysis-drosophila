# Transcription data analysis code from the [arXiv:1701.06079](https://arxiv.org/abs/1701.06079) manuscript.

To launch the code, open the `main.py` file. The lines of code are accompanied by comments for an easier understanding of the function. No installation is required to launch the code.

Every function located in a different file normally has a description at the top of the file. Most constants re-used throughout the analysis can be found and modified in the `constants.py`.

When running the `main.py` file from the terminal, comment out the first `try... except...` block, which provides code reloading for interactive environments like `Jupyter`.

## Input data file

To use the code, you need to supply the input data in a form of a table with pre-defined columns. The `.csv` file must be a comma-separated UTF-8 encoded file. An example of the input file called `example_data.csv` can be found in the `example` folder. The following columns are required to execute the code.

|Column | Description|
--- | ---
*time* | Time in seconds
*frame* | An integer corresponding to the frame number of the recorded video. Used to detect the time step
*intensity* | Fluorescence intensity of the transcription region in the nucleus, in a.u.
*ap* | Location of the region on the anterior-posterior axis of the embryo, where ap=0 corresponds to the embryo's head
*dataset* | Dataset name. Can be any string not containing slashes or backslashes
*dataset_id* | A unique sequential identifier with a one-to-one correspondence to the data set name. The counting starts from 0
*construct* | The name of the gene construct in the data set. Each data set must contain data from a single construct. In the data sets of the paper cited above, the construct takes one of 3 values: `bac`, `no_sh` or `no_pr`
*construct_id* | A unique sequential integer identifier for the constructs, starting at 0
*gene* | Gene name. Each data set must contain data from a single gene. In the data sets of the paper cited above, the construct takes one of 3 values: `hb` (*hunchback*), `sn` (*snail*) or `kn` (*knirps*)
*gene_id* | Same as `construct_id`, but for genes
*trace_id* | A sequential integer identifier for fluorescent traces (nuclei), starting at 0. Must be unique within each data set. Each data set may contain multiple traces

<!-- *nc* | Nuclear cycle, in which the current frame is recorded. An individual trace can span over multiple nuclear cycles. Only values in the range from 11 to 14 are processed in the code -->

## Requirements

Python 3

## Algorithm details

### Nuclear cycles detection

The nuclear cycles are detected in each data set by first selecting a background threshold, which separates the continuous trace into nuclear cycles above the threshold.
The default intensity threshold is defined by `default_intensity_threshold` in `constants.py`.
If a data set needs individual adjustments to the threshold, they must be put into the `intensity_thresholds` dictionary in `constants.py`.


The ncs are then numbered. By default, it is assumed that the last observed nc is nc14.
In some data sets, not all of the ncs have been recorded, or artefact ncs may be created by thresholding.
To correct for this problem, one may manually specify the number of the last nc in the data set by adding it into the `last_ncs` dictionary in `identify_ncs.py`.










