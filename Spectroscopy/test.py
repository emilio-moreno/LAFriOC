class Pal:
	def __init__(self):
		pass

class Friend(Pal):
	def __init__(self):
		pass


friend = Friend()
pal = Pal()
print(isinstance(pal, Friend))