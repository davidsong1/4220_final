#!/usr/bin/python
import os
import subprocess
with open("README.md", "w") as f:
    images = subprocess.run(["find", "images", "-type", "f"], stdout=subprocess.PIPE).stdout.decode('utf-8')
    f.write(images)
    with open("commands.txt", "r") as c:
        f.write(c.readline())
        f.write("yo")
        f.write(c.readline())
        c.close()
    f.close()
