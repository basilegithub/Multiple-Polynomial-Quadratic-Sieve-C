#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LINE 256

void trim(char *str)
{
    // remove trailing newline/whitespace
    char *end = str + strlen(str) - 1;
    while (end > str && (*end == '\n' || *end == ' ' || *end == '\t'))
        *end-- = '\0';
}

void parse_config(char* config_path, int* nb_cpu_sieve, int* flag_batch_smooth)
{
    FILE *file = fopen(config_path, "r");
    if (!file) {
        perror("Could not open config file");
    }

    char line[MAX_LINE];
    

    while (fgets(line, sizeof(line), file)) {
        trim(line);

        if (line[0] == '#' || strlen(line) == 0) // skip comments/empty lines
            continue;

        char *key = strtok(line, "=");
        char *value = strtok(NULL, "=");

        if (!key || !value) continue;

        if (strcmp(key, "nb_cpu_sieve") == 0) {
            *nb_cpu_sieve = atoi(value);
        } else if (strcmp(key, "flag_batch_smooth") == 0) {
            *flag_batch_smooth = atoi(value);
        }
    }

    fclose(file);
}
