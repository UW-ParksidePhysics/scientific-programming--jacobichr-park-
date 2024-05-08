def convert_units(value, from_units, to_units):
  """
  Converts a given value from one unit to another specified unit.

  Parameters:
  value (float): The numerical value to convert.
  from_units (str): The current units of the value ('cubic_bohr_per_atom', 'rydberg_per_atom', 'rydberg_per_cubic_bohr').
  to_units (str): The target units of the value ('cubic_angstroms_per_atom', 'electron_volts_per_atom', 'gigapascals').

  Returns:
  float: The value converted into the requested units.
  """
  # Define conversion factors
  conversions = {
      ('cubic_bohr_per_atom', 'cubic_angstroms_per_atom'): 0.14818471147216278,
      ('rydberg_per_atom', 'electron_volts_per_atom'): 13.605693122994,
      ('rydberg_per_cubic_bohr', 'gigapascals'): 14710.507848260711
  }

  # Check if the conversion is valid
  if (from_units, to_units) in conversions:
      # Perform the conversion
      converted_value = value * conversions[(from_units, to_units)]
      return converted_value
  else:
      raise ValueError(f"Conversion from {from_units} to {to_units} is not supported.")

if __name__ == "__main__":
  # Test cases
  print("Converting 1 cubic bohr per atom to cubic angstroms per atom:",
        convert_units(1, 'cubic_bohr_per_atom', 'cubic_angstroms_per_atom'))
  print("Converting 1 rydberg per atom to electron volts per atom:",
        convert_units(1, 'rydberg_per_atom', 'electron_volts_per_atom'))
  print("Converting 1 rydberg per cubic bohr to gigapascals:",
        convert_units(1, 'rydberg_per_cubic_bohr', 'gigapascals'))
