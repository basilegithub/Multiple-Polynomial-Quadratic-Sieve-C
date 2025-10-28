#ifndef CONFIG_H
#define CONFIG_H

void trim(char *str);
void parse_config(char* config_path, int* nb_cpu_sieve, int* flag_batch_smooth, int *flag_gaussian_elimination, int *flag_block_lanczos, size_t *block_size);

#endif // CONFIG_H