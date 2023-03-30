The workflow from raw data to pad vs. gap response diagrams relies on both automatic scripts and manual analysis.

An outline of the steps performed is as follows:

1. Demodulate raw data at the stimulus frequencies
   - raw data is in the `raw_data` directory
   - run `demodulation_script.m`
   - demodulated data is output into the `bw1000Hz` directory
2. Perform open/load impedance compensation
   - run `compensation_script.m`
   - compensated data is added to the same matfiles in the `bw1000Hz` directory
3. Automatic signal detection
   - run `new_find_padarr_signals.m`
   - cell transit events are detected and output as matfiles in the `detections` directory
4. Manual pruning of automatic detections
   - run `sort_padarr_signals.m`
   - specify which matfile to open
   - follow the prompts in the command window to classify each candidate detection as good or bad
   - good detections are saved as new matfiles in the `good_detections` directory
5. Fix misaligned detections
   - run `fix_good_detections`
   - follow the prompts in the command window to fix timestamps for starting edge and ending edge of each signal
   - fixed good detections are saved as new matfiles in the `fixed_detections` directory
6. Extract gap response and pad response values from signals
   - run first part of `analyze_fixed_detections`
   - this fits a deformable template to each signal
   - refined timestamps, fitting parameters, and filtering parameters are added to the good_detections struct array
   - run the last part of `analyze_fixed_detections`
   - this further refines the timestamps for the falling edge and rising edge of each pad response portion of the signal
   - the pad response and gap response values are extracted by averaging within a window between the time stamps
   - reponse values and further refined timestamps are added to the good_detections struct array
   - analyzed detections are saved as new matfiles in the `analyzed_detections` directory
7. Create scatter plots (pad vs. gap response diagrams)
   - run `make_scatter_plots.m`
   - select data using mask variable
   - performs fit of forward model to data to optimize CPE parameters (`q_vec` and `a_vec`)
   - plots fit lines over scatter plots
   - fitted parameters and saved as matfiles in `output` directory
   - plots are saved in `output` directory

