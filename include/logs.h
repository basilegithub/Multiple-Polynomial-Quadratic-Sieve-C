#ifndef LOGS_H
#define LOGS_H

void get_timestamp(char *buf, size_t size, const char *fmt);
void init_log(FILE **logfile);
void log_msg(FILE *logfile, const char *fmt, ...);
char *gmp_encode_string(const char *fmt, ...);
void log_gmp_msg(FILE *logfile, const char *fmt, ...);
void log_blank_line(FILE *logfile);

#endif // LOGS_H