CC=gcc
CFLAGS=-std=c11 -Wall -Werror -Wextra -g
SOURCES=s21_matrix.c
OBJECTS=$(SOURCES:.c=.o)
OS = $(shell uname)
TEST_FILE_NAME = s21_test_matrix.c

ifeq ($(OS), Darwin)
	TEST_LIBS=-lcheck
else
	TEST_LIBS=-lcheck -lsubunit -pthread -lrt -lm 
endif

all: s21_matrix.a

s21_matrix.a: clean $(OBJECTS)
	clang-format --style=google -i *.c *.h
	ar -rcs s21_matrix.a $(OBJECTS)
	rm -f *.o

test: $(TEST_FILE_NAME) s21_matrix.a
	$(CC) $(CFLAGS) $(TEST_FILE_NAME) $(SOURCES) -o test $(TEST_LIBS) -L. --coverage
	./test

gcov_report: test
	lcov -t "test" -o test.info -c -d .
	genhtml -o report test.info
	open report/index.html

clean:
	rm -rf *.o *.a *.so *.gcda *.gcno *.gch rep.info *.html *.css test report *.txt test.info test.dSYM proj

install_lcov:
	curl -fsSL https://rawgit.com/kube/42homebrew/master/install.sh | zsh
	brew install lcov

check: test
	cppcheck --enable=all --suppress=missingIncludeSystem --inconclusive --check-config *.c *.h
	rm -rf CPPLINT.cfg
	clang-format --style=google -n *.c *.h
	make test

leaks: test
	CK_FORK=no leaks --atExit --list -- ./test

ifeq ($(OS), Darwin)
	leaks --atExit -- test
else
	CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=RESULT_VALGRIND.txt ./test
endif
