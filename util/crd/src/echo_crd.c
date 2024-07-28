#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/crd.h"

struct rh1 h1;
struct rh2 h2;
struct rh3 h3;
struct rh4 h4;
struct rh5 h5;
struct rc0 c0;
struct rc1 c1;
struct rc2 c2;
struct rc3 c3;
struct rc4 c4;
struct rc5 c5;
struct rc6 c6;
struct rc7 c7;
struct rd10 d10;
struct rd11 d11;
struct rd12 d12;
struct rd20 d20;
struct rd21 d21;
struct rd30 d30;
struct rd40 d40;
struct rd40 d41;
struct rd42 d42;
struct rd50 d50;
struct rd60 d60;
struct rd00 d00;

int
main (argc, argv)
     int argc;
     char *argv[];
{
  char str[512];
  FILE *str_in;

  if (argc != 2)
    {
      printf ("Proper calling sequence is %s <input_file>\n", argv[0]);
      exit (0);
    }
  if ((str_in = fopen (argv[1], "r")) == NULL)
    {
      printf ("Could not open file %s\n", argv[1]);
      exit (0);
    }
  printf("in: %s\n", argv[1]);

  /* Copy and reformat data */
  while (fgets (str, 512, str_in) != NULL)
    {
      if (strncmp (str, "H1", 2) == 0 ||
          strncmp (str, "h1", 2) == 0)
        {
          read_h1 (str, &h1);
          printf ("\n--");
          printf ("\nRead H1 for type %s", h1.crd_literal);
        }
      else if (strncmp (str, "H2", 2) == 0 ||
               strncmp (str, "h2", 2) == 0)
        {
          read_h2 (str, &h2);
          printf ("\nRead H2 for station %s", h2.stn_name);
        }
      else if (strncmp (str, "H3", 2) == 0 ||
               strncmp (str, "h3", 2) == 0)
        {
          read_h3 (str, &h3);
          printf ("\nRead H3 for spacecraft %s", h3.target_name);
        }
      else if (strncmp (str, "H4", 2) == 0 ||
               strncmp (str, "h4", 2) == 0)
        {
          read_h4 (str, &h4);
          printf ("\nRead H4 for data type %i", h4.data_type);
        }
      else if (strncmp (str, "H5", 2) == 0 ||
               strncmp (str, "h5", 2) == 0)
        {
          read_h5 (str, &h5);
          printf ("\nRead H5 for prediction type %i", h5.prediction_type);
        }
      else if (strncmp (str, "H8", 2) == 0 ||
               strncmp (str, "h8", 2) == 0)
        {
          printf ("\nEnd of Session");
        }
      else if (strncmp (str, "H9", 2) == 0 ||
               strncmp (str, "h9", 2) == 0)
        {
          printf ("\n\nEnd of File");
        }
      else if (strncmp (str, "C0", 2) == 0 ||
               strncmp (str, "c0", 2) == 0)
        {
          read_c0 (str, &c0);
          printf ("\nRead C0 detail type %i", c0.detail_type);
        }
      else if (strncmp (str, "C1", 2) == 0 ||
               strncmp (str, "c1", 2) == 0)
        {
          read_c1 (str, &c1);
          printf ("\nRead C1 detail type %i", c1.detail_type);
        }
      else if (strncmp (str, "C2", 2) == 0 ||
               strncmp (str, "c2", 2) == 0)
        {
          read_c2 (str, &c2);
          printf ("\nRead C2 detail type %i", c2.detail_type);
        }
      else if (strncmp (str, "C3", 2) == 0 ||
               strncmp (str, "c3", 2) == 0)
        {
          read_c3 (str, &c3);
          printf ("\nRead C3 detail type %i", c3.detail_type);
        }
      else if (strncmp (str, "C4", 2) == 0 ||
               strncmp (str, "c4", 2) == 0)
        {
          read_c4 (str, &c4);
          printf ("\nRead C4 detail type %i", c4.detail_type);
        }
      else if (strncmp (str, "C5", 2) == 0 ||
               strncmp (str, "c5", 2) == 0)
        {
          read_c5 (str, &c5);
          printf ("\nRead C5 detail type %i", c5.detail_type);
        }
      else if (strncmp (str, "C6", 2) == 0 ||
               strncmp (str, "c6", 2) == 0)
        {
          read_c6 (str, &c6);
          printf ("\nRead C6 detail type %i", c6.detail_type);
        }
      else if (strncmp (str, "C7", 2) == 0 ||
               strncmp (str, "c7", 2) == 0)
        {
          read_c7 (str, &c7);
          printf ("\nRead C7 detail type %i", c7.detail_type);
        }
      else if (strncmp (str, "10", 2) == 0)
        {
          read_10 (str, &d10);
          printf ("\nRead range record %s", d10.sysconfig_id);
        }
      else if (strncmp (str, "11", 2) == 0)
        {
          read_11 (str, &d11);
          printf ("\nRead normal point record %s", d11.sysconfig_id);
        }
      else if (strncmp (str, "12", 2) == 0)
        {
          read_12 (str, &d12);
          printf ("\nRead range supplement record %s", d12.sysconfig_id);
        }
      else if (strncmp (str, "20", 2) == 0)
        {
          read_20 (str, &d20);
          printf ("\nRead meterological record");
        }
      else if (strncmp (str, "21", 2) == 0)
        {
          read_21 (str, &d21);
          printf ("\nRead meterological supplement record");
        }
      else if (strncmp (str, "30", 2) == 0)
        {
          read_30 (str, &d30);
          printf ("\nRead pointing angles record");
        }
      else if (strncmp (str, "40", 2) == 0)
        {
          read_40 (str, &d40);
          printf ("\nRead calibration record %s", d40.sysconfig_id);
        }
      else if (strncmp (str, "41", 2) == 0)
        {
          read_41 (str, &d41);
          printf ("\nRead calibration record II %s", d41.sysconfig_id);
        }
      else if (strncmp (str, "42", 2) == 0)
        {
          read_42 (str, &d42);
          printf ("\nRead calibration shot record %s", d42.sysconfig_id);
        }
      else if (strncmp (str, "50", 2) == 0)
        {
          read_50 (str, &d50);
          printf ("\nRead session statistics record %s", d50.sysconfig_id);
        }
      else if (strncmp (str, "60", 2) == 0)
        {
          read_60 (str, &d60);
          printf ("\nRead compatibility record %s", d60.sysconfig_id);
        }
      else if (strncmp (str, "00", 2) == 0)
        {
          read_00 (str, &d00);
          printf ("\nRead comment record %s", d00.comment);
        }
    }
  printf ("\n");
}
