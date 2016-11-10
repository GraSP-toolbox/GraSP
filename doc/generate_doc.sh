#!/bin/bash

CUR_DIR=$(pwd)
cd ..
DIRS=$(ls -d --indicator-style=none */ */*/)
cd $CUR_DIR

mkdir output_mwk output_github_md &> /dev/null
cp tuto.mwk output_mwk/
pandoc -f mediawiki -t markdown_github tuto.mwk -o output_github_md/tuto.md

for d in ./ $DIRS
do
    git ls-files ../$d --error-unmatch &> /dev/null
    if [[ $? -ne 0 ]]
    then
        continue
    fi
    TMP_FILES=$(ls ../${d}*m 2> /dev/null)
    if [[ $(echo $TMP_FILES | wc -w) -eq 0 ]]
    then
        continue;
    fi
    mkdir -p output_mwk/$d output_github_md/$d &> /dev/null
    FILE_LIST=""
    for f in $TMP_FILES
    do
        git ls-files $f --error-unmatch &> /dev/null
        if [[ $? -eq 0 ]]
        then
            FILE_LIST="$FILE_LIST $(echo $f | sed -e 's@.*/@@g')"
            cd ../$d
            OUTPUT_BASENAME=$(echo $f | sed -e 's@^\.\./@@')
            awk -f ${CUR_DIR}/matlabhelp2mediawiki.awk $(echo $f | sed -e 's@.*/@@g') > ${CUR_DIR}/output_mwk/${OUTPUT_BASENAME}.mwk
            pandoc -f mediawiki -t markdown_github ${CUR_DIR}/output_mwk/${OUTPUT_BASENAME}.mwk -o ${CUR_DIR}/output_github_md/${OUTPUT_BASENAME}.md
            cd ${CUR_DIR}
        fi
    done
    echo "== Files ==" > output_mwk/$d/index.mwk
    echo "" >> output_mwk/$d/index.mwk
    for f in $FILE_LIST
    do
        echo "* [[$(echo $f | sed -e 's/\.m$//')]]" >> output_mwk/$d/index.mwk
    done
    pandoc -f mediawiki -t markdown_github output_mwk/$d/index.mwk -o output_github_md/$d/index.md
done


for d in ./ $DIRS
do
    DIR_LIST=$(ls -d --indicator-style=none output_mwk/$d/*/ 2> /dev/null)
    if [[ $(echo $DIR_LIST | wc -w) -eq 0 ]]
    then
        continue;
    fi

    echo "" >> output_mwk/$d/index.mwk
    echo "== Directories ==" >> output_mwk/$d/index.mwk
    echo "" >> output_mwk/$d/index.mwk
    for sd in $DIR_LIST
    do
        echo "* [[$(echo $sd | sed -e 's@/$@@' | sed -e 's@.*/@@g')]]" >> output_mwk/$d/index.mwk
    done
    pandoc -f mediawiki -t markdown_github output_mwk/$d/index.mwk -o output_github_md/$d/index.md
done
