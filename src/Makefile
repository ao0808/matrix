CC = g++ -Wall -Werror -Wextra -std=c++17 

all: s21_matrix_oop.a

s21_matrix_oop.a:
	$(CC) -c s21_matrix_oop.cc -o s21_matrix_oop.o
	ar rcs s21_matrix_oop.a *.o
		ranlib s21_matrix_oop.a

test: s21_matrix_oop.a
	$(CC) test.cc -L. s21_matrix_oop.a -lcheck -lgtest -o test.out
	./test.out

clean:
	rm -rf *.o *.a *.out *.info report test.out.dSYM *.gcno

gcov_report: s21_matrix_oop.a
	$(CC) --coverage s21_matrix_oop.cc test.cc -lgtest s21_matrix_oop.a -L. s21_matrix_oop.a -o test.out
	./test.out
	lcov -t "my_test" -c -d ./ --output-file ./test.info
	genhtml -o report test.info
	rm -rf *.gcda *.gch *.gcno *.out *.o *.out
	open report/index.html

check: test
	leaks --atExit -- ./test.out

style:
	clang-format -i --style=Google *.cc *.h
	clang-format -n --style=Google *.cc *.h