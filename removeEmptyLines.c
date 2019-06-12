#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    const int maxline = 1024 * 1024;

    char* line = (char*)malloc(maxline);

    while (fgets(line, maxline, stdin)) {
        /* trim off newline if exists */
        if (line[strlen(line) - 1] == '\n' && strlen(line) == 1) {
            continue;
        }
        fputs(line, stdout);
    }

    return 0;
}
