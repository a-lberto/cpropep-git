# Bringing `cpropep` source to the current decade
## Disclaimer
All content is developed by the authors mentioned in the comments of the files.

I am keeping this copy for personal research purposes. Contact me if you are the author and wish for this to be removed.

```
/* Copyright (C) 2000                                                  */
/*    Antoine Lefebvre <antoine.lefebvre@polymtl.ca>                   */
/*    Mark Pinese <pinese@cyberwizards.com.au>                         */
/*                                                                     */
/* Licensed under the GPLv2                                            */
```

## Introduction
In my rocket propulsion studies I encountered a nifty tool called `cpropep` and have been fascinated by it ever since.

This fascination started with me developing a MatLab API to the compiled tool (which can be found at [thrust-team/mpropep](https://www.github.com/thrust-team/mpropep)) and ended up in the rabbit hole of me trying to reverse engineer the original source code and learn `bash`, `git` and ANSI C.

The following is a sort of log diary of the steps taken to create this repo.

## First try
At first I downloaded the the full repository (following the [instructions on the SourceForge rocketworkbench repository](http://rocketworkbench.cvs.sourceforge.net/) with the use of `cvs`), removed all CVS folders and tried to compile it loading the folders in MS Visual Studio. After some tricks and moving source and library files around the folders, the executable was working but the process was not scriptable and I did not learn much from it.

## Environment
To do these tasks I suggest a Unix environment. If you are on Windows I suggest WSL or a Virtual Machine.

To install WSL easily open an Admin Powershell `Win`+`X`+`A`, enable WSL and install for instance Debian:
```powershell
Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
Invoke-WebRequest -Uri https://aka.ms/wsl-debian-gnulinux -OutFile Debian.appx -UseBasicParsing
Add-AppxPackage .\Debian.appx 
```

Open the bash shell (open start menu and search "Debian") then install `git` and `git cvs`
```bash
sudo apt install git git-cvs
```

Make sure you have `curl` and `gpg`, add GitHub as a repository for apps, and then install `gh`
```bash
sudo apt install curl gpg
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | sudo gpg --dearmor -o /usr/share/keyrings/githubcli-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | sudo tee /etc/apt/sources.list.d/github-cli.list > /dev/null
sudo apt update
sudo apt install gh
```
Login to your GitHub account
```bash
gh auth login
```

## Unix way
I started looking for alternative ways to import a CVS repository into git, and found [a way](https://stackoverflow.com/a/11490134)

```bash
git cvsimport -C target-cvs -r cvs -k -vA authors-file.txt -d $CVSROOT module
```
which, if adapted to my case, becomes
```bash
mkdir rocketworkbench
cd rocketworkbench
git cvsimport -v -d :pserver:anonymous@a.cvs.sourceforge.net:/cvsroot/rocketworkbench $module -C $module
cd ..
```
where `$module`s are all the libraries needed to compile `cpropep`, listed in the `Makefile`:
- `cpropep`
- `libcpropep`
- `libnum`
- `libthermo`
- `libcompat`

Each folder of these modules has become a separate `git` repository, which are to be merged following [Eric Lee's advice](https://saintgimp.org/2013/01/22/merging-two-git-repositories-into-one-repository-without-losing-file-history/), although his `mv` commands did not work well for me.

I started by creating a new repo which would include all repos to be merged
```bash
mkdir cpropep-git
cd crpopep-git
git init
```
and by making a first commit with a dummy file
```bash
dir > deleteme.txt
git add .
git commit -m "Initial dummy commit"
```

Now, for every `$module`, a `git remote` will be added to the `cpropep-git` repo
```bash
for module in cpropep libnum libthermo libcpropep libcompat
do
    [...]
done
```

Inside the loop, I will do the following steps.

Add modules as remotes
```bash
git remote add -f $module ../rocketworkbench/$module
```

Merge the module history with the current repository [source](https://newbedev.com/how-to-automate-git-merge-to-not-prompt-to-confirm-the-commit-message)
```bash
git merge $module/master --allow-unrelated-histories --no-edit
```

Create subfolders to move the files to, move them and [delete the starting empty folders](https://unix.stackexchange.com/questions/8430/how-to-remove-all-empty-directories-in-a-subtree)
```bash
mkdir $module $module/doc/ $module/src/ $module/include/
mv doc/* $module/doc/
mv src/* $module/src/
mv include/* $module/include/
mv Makefile $module/Makefile
mv Makefile.win $module/Makefile.win
find . -type d -empty -delete
```

Commit the new module to the main repo:
```bash
git add .
git commit -m "Import $module files"
```

The `data` folder has a different structure and will be treated separately:
```bash
git remote add -f data ../rocketworkbench/data
git merge data/master --allow-unrelated-histories --no-edit
mkdir data
mv propellant.dat data/propellant.dat
mv thermo.dat data/thermo.dat
mv references.txt data/references.txt
git add .
git commit -m "Import data files"
```

Clean the dummy file inserted at the start
```bash
git rm ./deleteme.txt
git commit -m "Clean dummy file"
```

Create a GitHub repository and push the commit history [source](https://docs.github.com/en/github/importing-your-projects-to-github/importing-source-code-to-github/adding-an-existing-project-to-github-using-the-command-line)
```bash
gh repo create cpropep-git
git push --set-upstream origin master
```

## Result
The entire script is then:
```bash
mkdir rocketworkbench
cd rocketworkbench
for module in cpropep libnum libthermo libcpropep libcompat
do
    git cvsimport -v -d :pserver:anonymous@a.cvs.sourceforge.net:/cvsroot/rocketworkbench $module -C $module
done
cd ..

mkdir cpropep-git
cd cpropep-git
git init
dir > deleteme.txt
git add .
git commit -m "Initial dummy commit"

for module in cpropep libnum libthermo libcpropep libcompat
do
    git remote add -f $module ../rocketworkbench/$module
    git merge $module/master --allow-unrelated-histories --no-edit


    mkdir $module $module/doc/ $module/src/ $module/include/
    mv doc/* $module/doc/
    mv src/* $module/src/
    mv include/* $module/include/
    mv Makefile $module/Makefile
    mv Makefile.win $module/Makefile.win
    find . -type d -empty -delete

    git add .
    git commit -m "Import $module files"
done

git remote add -f data ../rocketworkbench/data
git merge data/master --allow-unrelated-histories --no-edit
mkdir data
mv propellant.dat data/propellant.dat
mv thermo.dat data/thermo.dat
mv references.txt data/references.txt
git add .
git commit -m "Import data files"

git rm ./deleteme.txt
git commit -m "Clean dummy file"

gh repo create cpropep-git
git push --set-upstream origin master
```
although running all commands together is not advised.

## The End
I did it! This repo is the result of merging the required libraries to compile `cpropep`. The last step is to provide continuous integration and let GitHub compile the source. (To be continued)
