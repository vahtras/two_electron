default:
	python setup.py build
install:
	python setup.py install
test:
	python -m pytest --cov=two --cov-report=html tests 
