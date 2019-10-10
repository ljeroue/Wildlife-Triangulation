# Wildlife-Triangulation

This set of functions estimate target location (potentially with error) using known
transmitter locations and bearings. Three estimation methods available for non-robust
and robust estimates following Lenth (1981).

Lenth, R.V. 1981.  On finding the source of a signal.
Technometrics 23:149-154. 

To run simply use the estTelemetry function.
estTelemetry is the main function which calls confEllipse and esestLocation

Note:
All data must be in the same coordinate reference system