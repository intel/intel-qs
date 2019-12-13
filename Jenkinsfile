pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
                echo 'Building..'
                echo '$USER'
                sh 'sudo usermod -a -G docker $USER'
                sh '$WORKSPACE/docker_build'
            }
        }
    }
}
