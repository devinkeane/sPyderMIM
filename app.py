from flask import Flask, request, render_template, send_from_directory, session, flash, redirect, url_for, jsonify
import os
import subprocess
import zipfile
import glob

app = Flask(__name__)
app.secret_key = 'super secret key'

# Get API key from environment variable
api_key = os.environ.get('OMIM_API_KEY')

@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        # Reset session variables
        session['download_ready'] = False
        session['download_file'] = ''

        mim_numbers = request.form['mim_numbers'].split('\n')
        project_name = request.form['project_name']

        with open('tmp.txt', 'w') as f:
            for mim_number in mim_numbers:
                f.write("%s\n" % mim_number)

        process = subprocess.run(['./sPyderMIM.sh', 'tmp.txt', api_key, project_name])

        # Get the most recently created folder that includes "SPYDERMIM_RESULTS" in its name
        output_folders = glob.glob('*SPYDERMIM_RESULTS*')
        output_folders.sort(key=os.path.getmtime, reverse=True)
        output_folder = output_folders[0]

        # Create a zip file from the output folder
        zipf = zipfile.ZipFile(f"{output_folder}.zip", 'w', zipfile.ZIP_DEFLATED)
        for root, dirs, files in os.walk(output_folder):
            for file in files:
                zipf.write(os.path.join(root, file),
                           os.path.relpath(os.path.join(root, file),
                                           os.path.join(output_folder, '..')))
        zipf.close()

        session['download_ready'] = True
        session['download_file'] = f"{output_folder}.zip"

        # Send a JSON response back to the AJAX request
        return jsonify(download_ready=session['download_ready'], download_file=session['download_file'])

    return render_template('home.html')

@app.route('/return-files/<filename>')
def return_files_tut(filename):
    session['download_ready'] = False
    return send_from_directory(directory='.', path=filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
