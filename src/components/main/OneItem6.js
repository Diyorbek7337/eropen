import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';


const OneItem6 = () => {
    const { t } = useTranslation()
    const [ block, setBlock ] = useState(false)
    return (
        <div className="block__item">
        <button className={block ? "block__title --active": "block__title"} type="button" data-spoller="" onClick={() => setBlock(!block)}>
        <i className="fas fa-first-aid"></i> {t('Home.releaseTitle')}
          <span className="spoller__arrow"></span>
        </button>
        <div className={block ? "block__text --active": "block__text"} hidden="">
        {t('Home.releaseText1')}
          <br/>
          {t('Home.releaseText2')}
        </div>
      </div>
    );
};

export default OneItem6;