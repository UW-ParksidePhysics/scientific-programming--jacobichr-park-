code_snippets = {
  'working': {
      'code': "numbers = {}\nnumbers[0] = -5\nnumbers[1] = 10.5",
      'explanation': "This snippet works because a dictionary in Python can have new keys assigned directly.",
      'fixed_code': None
  },
  'non_working': {
      'code': "other_numbers = []\nother_numbers[0] = -5\nother_numbers[1] = 10.5",
      'explanation': "This snippet doesn't work because a list's indices must be initialized before they can be assigned.",
      'fixed_code': "other_numbers = [None] * 2\nother_numbers[0] = -5\nother_numbers[1] = 10.5"
  }
}

def print_code_snippets(snippets):
  for snippet_name, details in snippets.items():
      print(f"Snippet: {snippet_name}")
      print("Code:")
      print(details['code'])
      print("Explanation:")
      print(details['explanation'])
      if details['fixed_code']:
          print("Fixed Code:")
          print(details['fixed_code'])
      print("-" * 50)

if __name__ == '__main__':
  print_code_snippets(code_snippets)
