# LWA_FRB
Scripts to process LWA beamformed data from CHIME triggers

- Voltage data in the form of raw files undergo coherent dedispersion and they are converted to fits files
- The fits files are read, data is read into the memory undergoes bandpass calibration, RFI flagging, incoherent dedispersion and candidate plots are made
- Slices of data around the expected FRB regions are saved into a .npz format
- Data slices for multiple events are stacked together based on the DM of the source

