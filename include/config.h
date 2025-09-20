#ifndef CONFIG_H
#define CONFIG_H

void trim(char *str);
void parse_config(char* config_path, int* nb_cpu_sieve, int* flag_batch_smooth);

#endif // CONFIG_H