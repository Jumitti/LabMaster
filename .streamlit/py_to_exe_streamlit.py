# MIT License
#
# Copyright (c) 2024 Minniti Julien
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import subprocess
import time

current_directory = os.path.dirname(os.path.abspath(__file__))
app_script_path = os.path.join(current_directory, "", "labmaster.py")
requirements_path = os.path.join(current_directory, "", "requirements.txt")

start_messages = [
    "------------------------------------------------------------------------------------------------------",
    "!                                        Welcome to LabMaster                                        !",
    "------------------------------------------------------------------------------------------------------",
    "LabMaster is designed to streamline and simplify the calculations essential for biological research. ",
    "Our software currently offers tools for transfection calculations, sample preparation for Western blot ",
    "analysis, reverse transcription calculations for RNA work, and Venn diagram generation.",
    "------------------------------------------------------------------------------------------------------",
    "Created by Minniti Julien - GitHub(https://github.com/Jumitti/LabMaster)",
    "MIT licence(https://github.com/Jumitti/LabMaster/blob/main/LICENSE)"
    ]
for message in start_messages:
    os.system(f"echo {message}")

time.sleep(2)

update_messages = "Installing/updating python prerequisites..."
os.system(f"echo {update_messages}")
subprocess.run(["pip", "install", "-r", requirements_path])

app_messages = "LabMaster Streamlit app running..."
subprocess.run(["streamlit", "run", app_script_path])
