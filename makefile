# hotplate make file
CC = gcc
CFLAGS = -pthread -O3
ver1 = v1_vanilla.c
ver2 = v2_fast_test.c
ver3 = v3_log_barrier.c
ver4 = v4_busy_wait.c

v1:
	$(CC) $(CFLAGS) $(ver1) -o hot
v2:
	$(CC) $(CFLAGS) $(ver2) -o hot
v3:
	$(CC) $(CFLAGS) $(ver3) -o hot
v4:
	$(CC) $(CFLAGS) $(ver4) -o hot
