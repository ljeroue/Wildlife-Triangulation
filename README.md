# Wildlife-Triangulation

description<< Estimates target location (potentially with error) using known
  # station locations and bearings and with user defined triangulation methods 
  # based on methods in: Lenth, R.V. 1981.  On finding the source of a signal.
  # Technometrics 23:149-154. Produces associated ArcGIS Shapefiles.

  # Returns:
  # A) A folder in user defined file location including shapefiles: 
  # 1) Triangulation_Stations,
  # 2) Triangulated_pts, 
  # 3) Bearings, and 
  # 4) Confidence_Region.
  # B) User can call a list composed of three elements: 
  # 1) A list contaning the von Mises estimate of location (vm), covariance 
  # matrix (qhat), weights (w), k, n, tuning constant, bearings, outcome, 
  # and method, 
  # 2) a data frame including location estimates and 
  # 3) a data frame including area of confidence region.

All data must be in the same coordinate reference system