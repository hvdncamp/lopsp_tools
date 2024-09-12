CC = gcc
TARGET_NAMES = apply_lopsp txt_to_lopsp decocode_to_lopsp lopspgen read_edgecode read_lopsp read_planarcode

all: $(TARGET_NAMES)

%: src/%.c graph_io.o lopsp_functions.o
	$(CC) -O3 -o $@ $^
	
%.o: src/lib/%.c
	$(CC) -O3 -o $@ -c $^

.PHONY: clean
clean:
	rm -rf obj/*.o $(TARGET_NAMES)
