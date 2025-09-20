#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <gmp.h>

// Helper: get current timestamp string
void get_timestamp(char *buf, size_t size, const char *fmt) {
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    strftime(buf, size, fmt, t);
}

// Initialize log file with run_YYYYMMDD_HHMMSS.txt
void init_log(FILE **logfile) {
    char filename[64];
    get_timestamp(filename, sizeof(filename), "./logs/run_%Y%m%d_%H%M%S.txt");

    *logfile = fopen(filename, "w");

    if (!*logfile) {
        perror("Could not open log file");
        exit(1);
    }
}

// Logging function (like printf)
void log_msg(FILE *logfile, const char *fmt, ...) {

    char timestamp[32];
    get_timestamp(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S");

    va_list args;

    // --- print to stdout ---
    va_start(args, fmt);
    printf("%s ", timestamp);
    vprintf(fmt, args);
    printf("\n");
    va_end(args);

    // --- write to logfile ---
    va_start(args, fmt);   // restart va_list before second use!
    fprintf(logfile, "%s ", timestamp);
    vfprintf(logfile, fmt, args);
    fprintf(logfile, "\n");
    fflush(logfile);
    va_end(args);

    va_end(args);
}

// GMP-aware logger
void log_gmp_msg(FILE *logfile, const char *fmt, ...) {
    if (!logfile) return;

    char timestamp[32];
    get_timestamp(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S");

    va_list args;

    // --- print to terminal ---
    va_start(args, fmt);
    printf("%s ", timestamp);
    gmp_vfprintf(stdout, fmt, args);
    printf("\n");
    va_end(args);

    // --- print to log file ---
    va_start(args, fmt);  // restart args!
    fprintf(logfile, "%s ", timestamp);
    gmp_vfprintf(logfile, fmt, args);
    fprintf(logfile, "\n");
    fflush(logfile);
    va_end(args);
}

// Print a completely blank line (no timestamp)
void log_blank_line(FILE *logfile) {
    if (!logfile) return;

    // print to stdout
    printf("\n");

    // print to logfile
    fprintf(logfile, "\n");
    fflush(logfile);
}