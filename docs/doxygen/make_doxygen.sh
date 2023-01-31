
#!/bin/bash

#root directory of project and config file for doxygen
parentname="ViennaEMC"
templatefile="doxygen_template.txt"
configfile="doxygen_config.txt"

#filename for html shortcut to Documentation outside of html folder
docmainfile="../index.html"

#EnsembleMC-directory is automatically found by searching current working directory
EnsembleMCdir=$(pwd | sed 's/\(.*ViennaEMC\).*/\1/' | sed -e 's/[\/&]/\\&/g')

echo "Using directory: ${EnsembleMCdir}"
#list of all commands and arguments that need to be changed for different computers
commands[0]="USE_MDFILE_AS_MAINPAGE"
commands[1]="INPUT"
commands[2]="EXCLUDE"
commands[3]="PROJECT_LOGO"
commands[4]="STRIP_FROM_PATH"
commands[5]="EXAMPLE_PATH"
commands[6]="USE_MATHJAX"

argument[0]="$EnsembleMCdir\/README.md"
argument[1]="$EnsembleMCdir" 
argument[2]="$EnsembleMCdir\/cmake $EnsembleMCdir\/tests $EnsembleMCdir\/helper $EnsembleMCdir\/docs $EnsembleMCdir\/dependencies" 
argument[3]="$EnsembleMCdir\/docs\/doxygen\/logo.png"
argument[4]="$EnsembleMCdir"
argument[5]="$EnsembleMCdir\/examples"
argument[6]="YES"

#create doxygen config file
cp $templatefile $configfile

#change the doxygen config file to the new parameters
for i in `seq 0 $(( ${#commands[*]} - 1 ))`; do
	command='sed -i "s/\(.*'${commands[$i]}' *=\).*/\1 '${argument[$i]}'/" '$configfile
	eval $command
done

#run doxygen to make the page
doxygen $configfile

#if shortcut outside the html folder does not exist, create it
if [ -e $docmainfile ]
then
	echo "$docmainfile already exists: Not creating."
else
echo '<!DOCTYPE html>
<html lang="en-US">
    <head>
        <meta charset="UTF-8">
        <meta http-equiv="refresh" content="1; url=doxygen/html/index.html">
        <script type="text/javascript">
            window.location.href = "doxygen/html/index.html"
        </script>
        <title>Page Redirection</title>
    </head>
    <body>
        If you are not redirected automatically, follow this <a href="doxygen/html/index.html">link to the manual</a>.
    </body>
</html>' > $docmainfile
echo "Created shortcut: $docmainfile"
fi