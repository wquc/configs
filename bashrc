alias ls='ls -v -Isnap --group-directories-first --color=auto --hide="*.pyc"'
alias ll='ls -alF'
alias la='ls -A'
alias l='ls -CF'
alias lll='ls -lht'
alias sl='ls'
alias e='exit'
alias h='history'
alias p='pwd'
alias ..='cd ..'
alias ...='cd ../..'
alias tf='tail -f'
alias get='readlink -f'
alias hg='history | grep'
alias vb='vi ~/.bashrc'
alias sb='source ~/.bashrc'
alias f='xdg-open'
alias CP='cp -rv'
alias RM='rm -rvf'
alias md='mkdir'
alias cdd='cd -P'
# requires [xclip] installed
alias pc="pwd | tee /dev/tty | tr -d '\n' | xclip -i"

# make backup for a file
bak() {
        cp $1{,.BAK}
}

# fetch a PDB file from RCSB site
fetch() {
       wget ""https://files.rcsb.org/download/$1.pdb""
}

# extract commonly used compressed files
extract () {
   if [ -f $1 ] ; then
       case $1 in
           *.tar.bz2)   tar xvjf $1    ;;
           *.tar.gz)    tar xvzf $1    ;;
           *.bz2)       bunzip2 $1     ;;
           *.rar)       unrar x $1     ;;
           *.gz)        gunzip $1      ;;
           *.tar)       tar xvf $1     ;;
           *.tbz2)      tar xvjf $1    ;;
           *.tgz)       tar xvzf $1    ;;
           *.zip)       unzip $1       ;;
           *.Z)         uncompress $1  ;;
           *.7z)        7z x $1        ;;
           *)           echo "don't know how to extract '$1'..." ;;
       esac
   else
       echo "'$1' is not a valid file!"
   fi
}

bind 'set completion-ignore-case on'
