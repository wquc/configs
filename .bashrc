LS_COLORS="$LS_COLORS:ow=103;30;01"                       # better color scheme for 777 files

HISTTIMEFORMAT=`echo -e "\033[1;30m" %F %T "\033[0m" `    # timestamp of history command with color

# Common
alias ls='ls -v --group-directories-first --color=auto --hide="*.pyc"'
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
alias md='mkdir'
alias cdd='cd -P'
alias CP='cp -rv'
alias RM='rm -rvf'
alias CHMOD='chmod -R'
alias pip='pip3'
alias python='python3'
alias jkl='jupyter lab' # need jupyterlab installed 

# For WSL
# alias o='explorer.exe'
# alias n='notepad.exe'
# export DISPLAY=127.0.0.1:0.0 

# make backup for a file
bak() {
    cp $1{,.BAK}
}

# make backup for a file with date (add %H%M for hours and minutes)
bakk() {
    cp $1{,.BAK-`date +%m%d%y`}
}

# make a directory by date
mdd() {
    mkdir BAK-`date +"%m-%d-%Y"`
}

# fetch a PDB file from RCSB site
fetch() {
    wget ""https://files.rcsb.org/download/$1.pdb""
}

# get full path of a file as input for rsync
rget() {
    host_addr=`hostname -I | cut -d ' ' -f1`
    file_path=`readlink -f $1`
    user_name=$USER
    echo ${user_name}@${host_addr}:${file_path}
}

# extract commonly used compressed files
x () {
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

# compress directory to tar.gz
z () {
    dirname=$1
    if [ -d $dirname ]; then
        tgzname=${dirname%%/}.tar.gz
        tar -czvf $tgzname $dirname
    fi
}

bind 'set completion-ignore-case on'
export PIPSRC="https://pypi.tuna.tsinghua.edu.cn/simple/"   # Domestic PIP source
# export PS1="\[\e[1;33m\]\u\[\e[m\]\[\e[1;33m\]@\[\e[m\]\[\e[1;33m\]\h\[\e[m\]:\W\\$ "
