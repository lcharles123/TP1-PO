prog: t1

t1:
	python3 main.py < teste0

todos:
	for i in `seq 1 4` ; do echo "teste e$$i" ; python3 main.py < "./tests/e$$i" || break ; done
	for i in `seq 1 7` ; do echo "teste $$i" ; python3 main.py < "./tests/$$i" || break ; done
