all: snapshot_integrity

snapshot_integrity: snapshot_integrity.cpp snapshot_integrity.h
	g++ -g -Wall -o snapshot_integrity snapshot_integrity.cpp snapshot_integrity.h -std=c++14

clean:
	rm -f snapshot_integrity