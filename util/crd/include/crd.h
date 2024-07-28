/*  Consolidated laser Ranging Data format (CRD)
 *       record and variable definitions for FORTRAN
 *       R. Ricklefs UT/CSR July 2007
 *  History:
 *  08/xx/07 - added H3 Target type (v0.27)
 *  11/26/07 - added H4 data quality alert
 *	       and #10 stop number
 *	       and #20 origin of values (v0.27) rlr.
 *  05/07/08 - Expand all configuration and data section character strings
 *             to allow up to 40 characters (plus NULL). 
 *             - Added detector channel to normalpoint (11) and calibration (40)
 *             records. 
 *             - Added field for 'crd' literal to H1.
 *             - Record '21' sky_clarity is not double rather than int.
 *               (v1.00 rlr)
 *  03/08/18 - Changes for CRD v2.00. rlr.
 *  06/26/19 - Added #42 calibration "shot" record for v2.00. rlr.
 *  07/xx/19 - Added NA_VALUE. v2.00. rlr.
 *  09/05/19 - Added NA_VALUEF. v2.01. rlr.
 *
 * ----------------------------------------------------------------------
*/

#ifndef _CRD_H
#define _CRD_H 1

/* Common name abbreviations:
 *  app = applied
 *  cdp = Crustal Dynamics Project (old NASA name)
 *  CofM = Center of Mass
 *  corr = correction or corrected
 *  est = estimated
 *  ind = indicator
 *  num = number
 *  osc = oscillator
 *  off = offset
 *  PmM = peak minus mean
 *  sic = (Goddard) Satellite ID Code
 *  stn = station
 *  SC = spacecraft
 *  sys = system
 *  utc = Universal Time Coordinated
 *  xcv = receive
 *  xmt = transmit
 */
#define NA_VALUE (-10000000)
#define NA_VALUEF (-1.0e30)

/* Ranging data header fields */
  /* H1 - format header */
  struct rh1 {
    char crd_literal[4];
    int format_version;
    int prod_year;
    int prod_mon;
    int prod_day;
    int prod_hour;
  };

  /* H2 - station header */
  struct rh2 {
    char stn_name[11];
    int cdp_pad_id;
    int cdp_sys_num;
    int cdp_occ_num;
    int stn_timescale;
    char stn_network[11];
  };

  /* H3 - spacecraft header */
  struct rh3 {
    char target_name[11];
    int ilrs_id;
    int sic;
    int norad;
    int SC_timescale;
    int target_type;
    int target_class;	// V2
    int target_loc;	// V2
  };

  /* H4 - Session header */
  struct rh4 {
    int data_type;
    int start_year;
    int start_mon;
    int start_day;
    int start_hour;
    int start_min;
    int start_sec;
    int end_year;
    int end_mon;
    int end_day;
    int end_hour;
    int end_min;
    int end_sec;
    int data_release;
    int refraction_app_ind;
    int CofM_app_ind;
    int xcv_amp_app_ind;
    int stn_sysdelay_app_ind;
    int SC_sysdelay_app_ind;
    int range_type_ind;
    int data_qual_alert_ind;
  };

  /* H5 - Prediction header */
  struct rh5 {			// V2
    int prediction_type;
    int year_of_century;
    char date_and_time[13];
    char prediction_provider[11];
    int  sequence_number;
  };

  /* Need indicators that these have been read? */
  /* H8 - End of Session footer */
  /* H9 - End of File footer */

/* Ranging data configuration fields (1 of n) */
  /* C0 - System Configuration Record */
  struct rc0 {
    int detail_type;
    double xmit_wavelength;
/**
    char sysconfig_id[5];
    char laserconfig_id[4];
    char detectorconfig_id[4];
    char timingconfig_id[4];
    char xponderconfig_id[4];
    char softwareconfig_id[4];
    char metconfig_id[4];
    char calconfig_id[4];
**/
    char config_ids[10][41];
  };

  /* C1 - Laser Configuration Record */
  struct rc1 {
    int detail_type;
    char laser_config_id[41];
    char laser_type[41];
    double prim_wavelength;	/* Primary wavelength of laser */
    double nom_fire_rate;	/* Nominal fire rate of laser */
    double pulse_energy;
    double pulse_width;
    double beam_div;
    int pulses_in_semitrain;	/* for multi-pulse systems */
  };

  /* C2 - Detector Configuration Record */
  struct rc2 {
    int detail_type;
    char detector_config_id[41];
    char detector_type[41];
    double app_wavelength;
    double qe;			/* quantum efficiency (in %) */
    double voltage;
    double dark_count;
    char output_pulse_type[41];
    double output_pulse_width;
    double spectral_filter;
    double spectral_filter_xmission;	/* % transmission of filter */
    double spatial_filter;
    char signal_proc[41];	/* signal processing algorithm or pgm name */
    double amp_gain;
    double amp_bandwidth;
    int amp_in_use;
  };

  /* C3 - Timing Configuration Record */
  struct rc3 {
    int detail_type;
    char timing_config_id[41];
    char time_source[41];
    char freq_source[41];
    char timer[41];
    char timer_serial_num[41];
    double epoch_delay_corr;
  };

  /* C4 - Transponder Configuration Record */
  struct rc4 {
    int detail_type;
    char xponder_config_id[41];
    long double est_stn_utc_offset;
    double est_stn_osc_drift;
    long double est_xponder_utc_offset;
    double est_xponder_osc_drift;
    long double xponder_clock_ref_time;
    int stn_off_drift_app_ind;
    int SC_off_drift_app_ind;
    int SC_time_simplified_ind;
  };

  /* C5 - Software Configuration Record */
  struct rc5 {			// V2
    int detail_type;
    char software_config_id[41];
    char tracking_software[41];
    char tracking_software_versions[41];
    char processing_software[41];
    char processing_software_versions[41];
  };

  /* C6 - Meteorology Configuration Record */
  struct rc6 {			// V2
    int detail_type;
    char met_config_id[41];
    char pressure_sensor_manufacturer[41];
    char pressure_sensor_model[41];
    char pressure_sensor_serial_num[41];
    char temperature_sensor_manufacturer[41];
    char temperature_sensor_model[41];
    char temperature_sensor_serial_num[41];
    char humidity_sensor_manufacturer[41];
    char humidity_sensor_model[41];
    char humidity_sensor_serial_num[41];
  };

  /* C7 - Calibration Target Configuration Record */
  struct rc7 {			// V2
    int detail_type;
    char calconfig_id[41];
    char target_name[41];
    double surveyed_target_dist;
    double survey_error;
    double other_fixed_delays;
    double pulse_energy;
    char processing_software[41];
    char processing_software_version[41];
  };

/* Ranging data fields */
  /* 10 - Range Record */
  struct rd10 {
    long double sec_of_day;
    long double time_of_flight;
    char sysconfig_id[41];
    int epoch_event;
    int filter_flag;
    int detector_channel;
    int stop_number;
    int xcv_amp;
    int xmt_amp;	// V2
  };

  /* 11 - Normal Point Record */
  struct rd11 {
    long double sec_of_day;
    long double time_of_flight;
    char sysconfig_id[41];
    int epoch_event;
    double np_window_length;
    int num_ranges;
    double bin_rms;
    double bin_skew;
    double bin_kurtosis;
    double bin_PmM;
    double return_rate;
    int detector_channel;
    double signal_to_noise;	// V2
  };

  /* 12 - Range Supplement Record */
  struct rd12 {
    long double sec_of_day;
    char sysconfig_id[41];
    double refraction_corr;
    double target_CofM_corr;
    double nd_value;
    double time_bias;
    double range_rate;	// V2
  };

  /* 20 - Meteorological Record */
  struct rd20 {
    long double sec_of_day;
    double pressure;
    double temperature;
    double humidity;
    int value_origin;
  };

  /* 21 - Meteorological Supplement Record */
  struct rd21 {
    long double sec_of_day;
    double wind_speed;
    double wind_direction;
    char weather_conditions[41]; // V2
    int visibility;
    double sky_clarity;
    int atmospheric_seeing;
    int cloud_cover;
    double sky_temperature;	// V2
  };

  /* 30 - Pointing Angles Record */
  struct rd30 {
    long double sec_of_day;
    double azimuth;
    double elevation;
    int direction_ind;
    int angle_origin_ind;
    int refraction_corr_ind;
    double azimuth_rate;	// v2
    double elevation_rate;	// V2
  };

  /* 40 - Calibration Record */
  struct rd40 {
    long double sec_of_day;
    int type_of_data;
    char sysconfig_id[41];
    int num_points_recorded;
    int num_points_used;
    double one_way_target_dist;
    double cal_sys_delay;
    double cal_delay_shift;
    double cal_rms;
    double cal_skew;
    double cal_kurtosis;
    double cal_PmM;
    int cal_type_ind;
    int cal_shift_type_ind;
    int detector_channel;
    int cal_span;
    double cal_return_rate;
  };

  /* 42 - Calibration "Shot" Record */
  struct rd42 {	// V2
    long double sec_of_day;
    long double time_of_flight;
    char sysconfig_id[41];
    char calconfig_id[41];
    double other_variable_delays;
    int type_of_data;
    int cal_type_ind;
    int filter_flag;
    int detector_channel;
    int stop_number;
    int cal_span;
    int xcv_amp;
    int xmt_amp;
  };

  /* 50 - Session Statistics Record */
  struct rd50 {
    char sysconfig_id[41];
    double sess_rms;
    double sess_skew;
    double sess_kurtosis;
    double sess_PmM;
    int data_qual_ind;
  };

  /* 60 - Compatibility Record */
  // OBSOLETE -- V2
  struct rd60 {
    char sysconfig_id[41];
    int sys_change_ind;
    int sys_config_ind;
  };

  /* 9X - User Defined Record */
  struct rd9x {
    /********** 
        Add userdefined record types and fields here 
    **********/
  };

  /* 00 - Comment Record */
  struct rd00 {
    char comment[81];
  };

int read_h1 (char *, struct rh1 *);
int read_h2 (char *, struct rh2 *);
int read_h3 (char *, struct rh3 *);
int read_h4 (char *, struct rh4 *);
int read_h5 (char *, struct rh5 *);
int read_h8 (char*);
int read_h9 (char*);
int read_c0 (char *, struct rc0 *);
int read_c1 (char *, struct rc1 *);
int read_c2 (char *, struct rc2 *);
int read_c3 (char *, struct rc3 *);
int read_c4 (char *, struct rc4 *);
int read_c5 (char *, struct rc5 *);
int read_c6 (char *, struct rc6 *);
int read_c7 (char *, struct rc7 *);
int read_10 (char *, struct rd10 *);
int read_11 (char *, struct rd11 *);
int read_12 (char *, struct rd12 *);
int read_20 (char *, struct rd20 *);
int read_21 (char *, struct rd21 *);
int read_30 (char *, struct rd30 *);
int read_40 (char *, struct rd40 *);
int read_41 (char *, struct rd40 *);
int read_42 (char *, struct rd42 *);
int read_50 (char *, struct rd50 *);
int read_60 (char *, struct rd60 *);
int read_9x (char *, struct rd9x *);
int read_00 (char *, struct rd00 *);

int write_h1 (FILE *, struct rh1);
int write_h2 (FILE *, struct rh2);
int write_h3 (FILE *, struct rh3);
int write_h4 (FILE *, struct rh4);
int write_h5 (FILE *, struct rh5);
int write_h8 (FILE *);
int write_h9 (FILE *);
int write_c0 (FILE *, struct rc0);
int write_c1 (FILE *, struct rc1);
int write_c2 (FILE *, struct rc2);
int write_c3 (FILE *, struct rc3);
int write_c4 (FILE *, struct rc4);
int write_c5 (FILE *, struct rc5);
int write_c6 (FILE *, struct rc6);
int write_c7 (FILE *, struct rc7);
int write_10 (FILE *, struct rd10);
int write_11 (FILE *, struct rd11);
int write_12 (FILE *, struct rd12);
int write_20 (FILE *, struct rd20);
int write_21 (FILE *, struct rd21);
int write_30 (FILE *, struct rd30);
int write_40 (FILE *, struct rd40);
int write_41 (FILE *, struct rd40);
int write_42 (FILE *, struct rd42);
int write_50 (FILE *, struct rd50);
int write_60 (FILE *, struct rd60);
int write_9x (FILE *, struct rd9x);
int write_00 (FILE *, struct rd00);


#endif /* crd.h  */


