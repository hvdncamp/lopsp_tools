CC = gcc
TARGET_NAMES = apply_lopsp txt_to_lopsp decocode_to_lopsp lopspgen read_edgecode read_lopsp

all: $(TARGET_NAMES)

%: src/%.c graph_io.o lopsp_functions.o
	$(CC) -o $@ $^
	
%.o: src/lib/%.c
	$(CC) -o $@ -c $^

.PHONY: clean
clean:
	rm -rf obj/*.o $(TARGET_NAMES)
