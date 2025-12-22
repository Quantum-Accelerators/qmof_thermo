from logging import basicConfig, getLevelName

basicConfig(filename=_settings.LOG_FILENAME, level=getLevelName(_settings.LOG_LEVEL))
