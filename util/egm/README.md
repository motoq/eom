egm
===

Earth Gravity Model Utility
---------------------------

This utility application reads in an earth gravity model file from:

[Office of Geomatics](https://earth-info.nga.mil)

The file is parsed to the desired degree and order.  The coefficients
are unnormalized and sorted, grouping equal order terms together.  If an
output file is specified, the terms are written to a C++ compatible
header file with separate, equally sized arrays, containing degree,
order, cosine, and sine terms.
