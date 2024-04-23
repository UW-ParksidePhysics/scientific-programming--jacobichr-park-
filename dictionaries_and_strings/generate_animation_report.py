import matplotlib.pyplot as plt
import numpy as np


def parse_animation_code(code_filename):
  with open(code_filename, 'r') as file:
      code_content = file.read()
  preformatted_code = '<pre>' + code_content + '</pre>'
  return preformatted_code


def format_section_header(header_string):
  return '<h1>' + header_string + '</h1>'


def write_html_report(report_string, report_filename):
  with open(report_filename, 'w') as file:
      file.write(report_string)

if __name__ == '__main__':
  # Read and format the code
  animation_code_html = parse_animation_code('dictionaries_and_strings/orbit_gif_creating_code.py')

  # Generate plots
  times = [np.pi/4, np.pi/2, 3*np.pi/4]  # Example times
  a, b, omega = 1, 1, 1  # Same constants as in the orbit animation code

  fig, axs = plt.subplots(1, 3, figsize=(15, 5))
  for i, time in enumerate(times):
      x = a * np.cos(omega * time)
      y = b * np.sin(omega * time)
      axs[i].plot(x, y, 'bo')
      axs[i].set_xlim(-2, 2)
      axs[i].set_ylim(-2, 2)
      axs[i].set_title(f't = {time:.2f}')
      axs[i].set_aspect('equal')
  plt.savefig('plots.png')

  # Prepare the HTML content
  html_content = """
  <html>
  <head><title>Orbit Animation Report</title></head>
  <body>
  """ + format_section_header('Orbit Simulation Code') + """
  """ + animation_code_html + """
  """ + format_section_header('Orbit at Three Different Times') + """
  <img src='plots.png' alt='Orbit Plots'>
  """ + format_section_header('Orbit Animation GIF') + """
  <img src='orbit_animation.gif' alt='Orbit Animation'>
  </body>
  </html>
  """

  # Write the HTML report
  write_html_report(html_content, 'animation_report.html')
