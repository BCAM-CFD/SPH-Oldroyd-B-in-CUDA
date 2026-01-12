// Following two lines avoid to process the content of the file more than once
#ifndef CONFIG_H
#define CONFIG_H

// Comment or uncomment to decide if real means float or double
// Type and format for double
/* #define real double */
/* #define REAL_FMT "%.8lf" */
// Type and format for real
#define real float
#define REAL_FMT "%.8f"

#endif // CONFIG_H
#define MAX_COLLOIDS 1000
#define MAX_VALUES_PER_COLLOID 3
