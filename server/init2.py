# -*- coding: utf-8 -*-

# from werkzeug.contrib.profiler import ProfilerMiddleware
from flask import Flask, request, g, render_template



app = Flask(__name__, static_path='/static')	
app.config['DEBUG'] = True


@app.route('/survey')
def form():
	return render_template('survey.html')


# Execute the main program
if __name__ == '__main__':
	app.run(host='0.0.0.0',port=1905, debug=True)
